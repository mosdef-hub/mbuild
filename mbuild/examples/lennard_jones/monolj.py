import mbuild as mb
from mbuild import clone


class MonoLJ(mb.Compound):
    def __init__(self):
        super(MonoLJ, self).__init__()
        lj_proto = mb.Particle(name='LJ', pos=[0, 0, 0])
        
        pattern = mb.Grid3DPattern(5, 5, 5)
        pattern.scale(5)
        
        for pos in pattern:
            lj_particle = clone(lj_proto)
            mb.translate(lj_particle, pos)
            self.add(lj_particle)
            print(pos)


if __name__ == '__main__':
    monolj = MonoLJ()
    print(monolj)
