__author__ = 'sallai'

class Prototype(object):
    prototypes = dict()

    def __init__(self, kind, **kwargs):
        super(Prototype, self).__init__()

        if not Prototype.prototypes.has_key(kind):
             Prototype.prototypes[kind] = dict()
        else:
            print "warning: overriding prototype for "+kind

        for key, value in kwargs.iteritems():
            Prototype.prototypes[kind][key] = value


        # for key, value in kwargs.iteritems():
        #     setattr(self, key, value)
        #
        # if not Prototype.prototypes.has_key(kind):
        #     Prototype.prototypes[kind] = self
        # else:
        #     raise Exception("Prototype already defined for "+kind)

    # @staticmethod
    # def get(kind):
    #     return Prototype.prototypes[kind]

    @staticmethod
    def getAttr(kind, key, default=None):
        if not Prototype.prototypes.has_key(kind):
            return default

        if not Prototype.prototypes[kind].has_key(key):
            return default

        return Prototype.prototypes[kind][key]

# define atom kind prototypes containing radius and color
Prototype('C', radius=1.7, color='teal')
Prototype('C2', radius=1.7, color='teal')
Prototype('CB', radius=1.7, color='teal')
Prototype('H', radius=1.2, color='white')
Prototype('HB', radius=1.2, color='white')
Prototype('O', radius=1.52, color='red')
Prototype('O1', radius=1.52, color='red')
Prototype('F', radius=1.35, color='pink')
Prototype('N', radius=1.55, color='blue')
Prototype('P', radius=1.8, color='orange')
Prototype('Si', radius=2.10, color='yellow')
Prototype('Si1', radius=2.10, color='yellow')

