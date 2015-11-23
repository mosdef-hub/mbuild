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

        self._com = None

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

    @property
    def root_compound(self):
        """Return the root of the compound hierarchy. """
        for root_of_hierarchy in self.ancestors():
            pass
        return root_of_hierarchy

    @property
    def tier_names(self):
        """Return the names of the tiers that this compound is tagged with. """
        root = self.root_compound
        if root.tiers is None:
            raise ValueError('No tiers have been tagged in the root of this'
                             ' hierarchy, {}.'.format(root))

        all_names = list()
        for tier_name, parts in root.tiers.items():
            if self in parts:
                all_names.append(tier_name)
        return all_names

    def _clone(self, clone_of=None, root_container=None):
        raise NotImplementedError
