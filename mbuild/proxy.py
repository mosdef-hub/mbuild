from collections import OrderedDict

from mbuild.compound import Compound


class Proxy(Compound):

    def __init__(self, compound):
        if compound.name == 'G':
            name = 'G'
        else:
            name = compound.name + ' (proxy) '
        super(Proxy, self).__init__(name=name)

        self.wrapped = compound
        self.children = None
        self.labels = None
        self.parent = None
        self.referrers = set()
        self.index = None
        self.bond_graph = None

    def proxy_for(self):
        if hasattr(self.wrapped, 'wrapped'):
            return self.wrapped.proxy_for()
        else:
            return self.wrapped.__class__

    @property
    def pos(self):
        return self.wrapped.pos

    @pos.setter
    def pos(self, value):
        self.wrapped.pos = value

    def __getattr__(self, attr):
        return getattr(self.wrapped, attr)


def is_leaf(what):
    return hasattr(what, 'parts') and not what.children


def create_proxy(real_thing, memo=None, particle_classes=None):
    if memo is None:
        memo = OrderedDict()

    if particle_classes is None:
        particle_classes = []

    proxy = _create_proxy_compounds(real_thing, memo, particle_classes)
    _create_proxy_bonds(real_thing, memo, particle_classes)
    _create_proxy_labels(real_thing, memo)

    return proxy


def _create_proxy_compounds(real_thing, memo, particle_classes):
    proxy = Proxy(real_thing)
    memo[real_thing] = proxy

    if not type(real_thing) in particle_classes:
        if not is_leaf(real_thing): # recurse only if it has parts
            # recursively create proxies for parts (we'll do labels later)
            for part in real_thing.children:
                part_proxy = _create_proxy_compounds(part, memo, particle_classes)
                proxy.add(part_proxy)

    return proxy


def _proxy_of(real_thing, memo):
    if real_thing in memo:
        return memo[real_thing]
    else:
        return _proxy_of(real_thing.parent, memo)


def _create_proxy_bonds(real_thing, memo, leaf_classes):
    proxy = memo[real_thing]

    if type(real_thing) in leaf_classes or is_leaf(real_thing):
        # it is a leaf of the proxy, so we don't recurse
        pass
    else:
        # recurse
        for part in real_thing.children:
            _create_proxy_bonds(part, memo, leaf_classes)

    # check if there's a contained bond that needs to be added to the proxy
    if hasattr(real_thing, 'contained_bonds'):
        for a1, a2 in real_thing.bonds:
            pa1 = _proxy_of(a1, memo)
            pa2 = _proxy_of(a2, memo)
            if pa1 != pa2:  # do not add internal bonds
                proxy.add_bond((pa1, pa2))


def _create_proxy_labels(real_thing, memo):
    if not is_leaf(real_thing):
        # create labels
        for label, part in real_thing.labels.items():
            if isinstance(part, list):
                # TODO support lists with labels
                continue
            if part in memo:
                memo[real_thing].labels[label] = memo[part]

        # recurse
        for part in real_thing.children:
            _create_proxy_labels(part, memo)


if __name__ == '__main__':
    from mbuild.examples import Ethane
    c = Ethane()
    # c = create_proxy(c)
    # p = create_proxy(c, leaf_classes=mbuild.lib.moieties.ch3.CH3)
    p = create_proxy(c)

    print("Top level proxy object: {}".format(p))

    print("Parts of top level proxy object:")
    for part in p.children:
        print(" {}".format(part))

    print("Leaves of top level proxy object:")
    for particle in p.particles():
        print(" {}".format(particle.name))
