class Part(object):
    """A base class that designates a part in a hierarchy. ""'

    Attributes
    ----------
    parent : mb.Compound
        The parent Compound that contains this part. Can be None if this
        compound is the root of the containment hierarchy.
    referrers : set
        Other compounds that reference this part with labels.
    tier : str
        Name of the user defined tier that this part belongs too.
    com : np.ndarray, shape=(3,)
        The center of mass of this part.
    """
    def __init__(self):
        super(Part, self).__init__()
        self.parent = None
        # The set of other compounds that reference this compound with labels.
        self.referrers = set()

        self._tier = None
        self._com = None

    @property
    def tier(self):
        return self._tier

    @tier.setter
    def tier(self, tier_name):
        self._tier = tier_name

    @property
    def com(self):
        return self._com

    @com.setter
    def com(self, center_of_mass):
        self._com = center_of_mass

    def ancestors(self):
        """Generate all ancestors of the Compound recursively. """
        yield self.parent
        if self.parent is not None:
            for ancestor in self.parent.ancestors():
                yield ancestor

    def _clone(self, clone_of=None, root_container=None):
        raise NotImplementedError
