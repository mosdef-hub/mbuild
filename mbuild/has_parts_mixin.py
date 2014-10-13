from collections import OrderedDict
from mbuild.mbase import MBase
from mbuild.orderedset import OrderedSet
from mbuild.part_mixin import PartMixin

__author__ = 'sallai'


class HasPartsMixin(object):
    def __init__(self, *args, **kwargs):
        super(HasPartsMixin, self).__init__()

        # Contains all child parts. Parts can be Atom, Bond or Compound.
        self.parts = OrderedSet()

        # Labels to Compound/Atom mappings. These do not necessarily need not
        # be in self.parts.
        self.labels = OrderedDict()

    def _yield_parts(self, part_type):
        """Generate all parts of a specified type in the Compound recursively.

        Args:
            part_type (Atom, Bond, Compound): The type of parts to yield.

        Yields:
            part/subpart (Atom, Bond, Compound): A part in the hierarchy
                matching the specified type.
        """
        for part in self.parts:
            # Parts local to the current Compound.
            if isinstance(part, part_type):
                yield part
            # Parts in further down the hierarchy.
            if isinstance(part, HasPartsMixin):
                for subpart in part._yield_parts(part_type):
                    yield subpart

    def add(self, new_part, label=None, containment=True, replace=False):
        """Add a part to the Compound.

        Note:
            This does not necessarily add the part to self.parts but may
            instead be used to add a reference to the part to self.labels. See
            'containment' argument.

        Args:
            new_part (Atom, Bond or Compound): The object to be added to this
                Compound.
            label (str, optional): A descriptive string for the part.
            containment (bool, optional):
            replace (bool, optional):
            inherit_periodicity (bool, optional):
        """
        assert isinstance(new_part, (PartMixin, list, tuple, set))
        if containment:
            # Support batch add via lists, tuples and sets.
            if isinstance(new_part, (list, tuple, set)):
                for elem in new_part:
                    assert(elem.parent is None)
                    self.add(elem)
                    elem.parent = self
                return
            # Is the following assertion necessary? Seems to be handled above.
            assert(new_part.parent is None)
            self.parts.add(new_part)
            new_part.parent = self

        # Add new_part to labels. Does not currently support batch add.
        assert isinstance(new_part, PartMixin)

        if not containment and label is None:
            label = '_{0}[$]'.format(new_part.__class__.__name__)

        if label is not None:
            if label.endswith("[$]"):
                label = label[:-3]
                if not label in self.labels:
                    self.labels[label] = []

                label_pattern = label + "[{}]"

                count = len(self.labels[label])
                self.labels[label].append(new_part)
                label = label_pattern.format(count)

            if not replace and label in self.labels:
                raise Exception("Label {0} already exists in {1}".format(label, self))
            else:
                self.labels[label] = new_part

        new_part.referrers.add(self)


    def remove(self, objs_to_remove):
        """Remove a part (Atom, Bond or Compound) from the Compound by value.

        Args:
            objs_to_remove (set of parts): All objects to be removed from the
                hierarchy. If this is not a set, it will be cast to one to
                remove duplicates.
        """
        if not isinstance(objs_to_remove, (list, tuple, set)):
            objs_to_remove = [objs_to_remove]

        objs_to_remove = set(objs_to_remove)

        if len(objs_to_remove) == 0:
            return

        intersection = objs_to_remove.intersection(self.parts)
        self.parts.difference_update(intersection)
        objs_to_remove.difference_update(intersection)

        for removed_part in intersection:
            self.post_remove(removed_part)

        # Remove the part recursively from sub-components.
        for part in self.parts:
            if isinstance(part, HasPartsMixin) and len(objs_to_remove) > 0:
                part.remove(objs_to_remove)


    def post_remove(self, removed_part):
        removed_part.parent = None
        # Remove labels in the hierarchy pointing to this part.
        referrers_to_remove = set()
        for referrer in removed_part.referrers:
            if not removed_part in referrer.ancestors():
                for label, referred_part in referrer.labels.items():
                    if referred_part is removed_part:
                        del referrer.labels[label]
                        # TODO:
                        # if label matches some_label[num], check if label some_label is an empty array, and remove it
                        referrers_to_remove.add(referrer)
        removed_part.referrers.difference_update(referrers_to_remove)

        # Remove labels in this part pointing into the hierarchy.
        labels_to_delete = []
        if isinstance(removed_part, HasPartsMixin):
            for label, part in removed_part.labels.items():
                if not removed_part in part.ancestors():
                    part.referrers.remove(removed_part)
                    labels_to_delete.append(label)
        for label in labels_to_delete:
            del removed_part.labels[label]

        # # If removing an atom, make sure to remove the bonds it's part of.
        # if isinstance(removed_part, Atom):
        #     for bond in removed_part.bonds:
        #         if bond.parent is not None:
        #             bond.parent.remove(bond)

    def __getattr__(self, attr):
        assert "labels" != attr, "HasPartsMixin __init__ never called. Make sure to call super().__init__() in the __init__ method of your class."
        if attr in self.labels:
            return self.labels[attr]
        else:
            raise AttributeError("'{}' object has no attribute '{}'".format(self.__class__.__name__, attr))
