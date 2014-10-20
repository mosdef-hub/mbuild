__author__ = 'sallai'


class PartMixin(object):
    """A base class that designates a part in a hierarchy. ""'

    Attribute:
        parent (Compound): The parent Compound that contains this part Can be
            None if this compound is the root of the containment hierarchy.
    """
    def __init__(self, *args, **kwargs):
        super(PartMixin, self).__init__()

        self.parent = None

        # The set of other compounds that reference this compound with labels.
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
