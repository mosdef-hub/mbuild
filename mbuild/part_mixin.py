from mbuild.mbase import MBase

__author__ = 'sallai'

class PartMixin(object):
    def __init__(self, *args, **kwargs):
        # The parent compound that contains this compound (can be None if this compound is the root of the
        # containment hierarchy.)

        super(PartMixin, self).__init__()

        self.parent = None

        # The set of other compounds that reference this compound with labels
        self.referrers = set()

    def ancestors(self):
        """Generate all ancestors of the Compound recursively.

        Yields:
            ancestor (Compound): A Compound one or more levels higher in the
                hierarchy.

        """
        yield self.parent
        if self.parent is not None:
            for ancestor in self.parent.ancestors():
                yield ancestor
