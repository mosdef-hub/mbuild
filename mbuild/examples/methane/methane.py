import mbuild as mb


class Methane(mb.Compound):
    def __init__(self):
        super(Methane, self).__init__()
        carbon = mb.Particle(name='C')
        self.add(carbon, label='C')

        hydrogen = mb.Particle(name='H', pos=[0.1, 0, -0.07])
        self.add(hydrogen, label='HC[$]')

        self.add_bond((self.C, self.HC[0]))

        self.add(mb.Particle(name='H', pos=[-0.1, 0, -0.07]), label='HC[$]')
        self.add(mb.Particle(name='H', pos=[0, 0.1, 0.07]), label='HC[$]')
        self.add(mb.Particle(name='H', pos=[0, -0.1, 0.07]), label='HC[$]')

        self.add_bond((self.C, self.HC[1]))
        self.add_bond((self.C, self.HC[2]))
        self.add_bond((self.C, self.HC[3]))


def main():
    methane = Methane()
    return methane

if __name__ == "__main__":
    methane = main()
    import ipdb; ipdb.set_trace()
    methane.visualize()
