import mbuild as mb
from mbuild import clone


class MonoLJ(mb.Compound):
    def __init__(self):
        super(MonoLJ, self).__init__()
        lj_particle1 = mb.Particle(name='C')
        self.add(lj_particle1)

        lj_particle2 = clone(lj_particle1)
        pos = [1, 0, 0]
        mb.translate(lj_particle2, pos)
        self.add(lj_particle2)


if __name__ == '__main__':
    monolj = MonoLJ()
    print(monolj)
    monolj.save('mono.mol2')
    monolj.visualize()
