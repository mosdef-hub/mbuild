__author__ = 'sallai'


class PartMixin(object):
    """A base class that designates a part in a hierarchy. ""'

    Attributes
    ----------
    parent : mb.Compound
        The parent Compound that contains this part. Can be None if this
        compound is the root of the containment hierarchy.
    referrers : set
        Other compounds that reference this part with labels.

    """
    def __init__(self):
        super(PartMixin, self).__init__()
        self.parent = None
        # The set of other compounds that reference this compound with labels.
        self.referrers = set()

    def ancestors(self):
        """Generate all ancestors of the Compound recursively. """
        yield self.parent
        if self.parent is not None:
            for ancestor in self.parent.ancestors():
                yield ancestor
