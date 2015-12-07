import mbuild as mb
from mbuild import clone


class MonoLJ(mb.Compound):
    def __init__(self):
        super(MonoLJ, self).__init__()
        lj_particle1 = mb.Particle(name='LJ')
        self.add(lj_particle1)

        lj_particle2 = clone(lj_particle1)
        pos = [1, 0, 0]
        mb.translate(lj_particle2, pos)
        self.add(lj_particle2)


if __name__ == '__main__':
    monolj = MonoLJ()

    colors = {'LJ': {'color': 0xbfbfbf, 'radius': 1.6}}
    monolj.visualize(element_properties=colors)
