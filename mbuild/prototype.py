from _warnings import warn

__author__ = 'sallai'

# TODO: Incorporate into plotting functions.
class Prototype(object):
    prototypes = dict() # kind:attribute_dict

    def __init__(self, kind, **kwargs):
        super(Prototype, self).__init__()

        if not Prototype.prototypes.has_key(kind):
             Prototype.prototypes[kind] = dict()

        for key, value in kwargs.iteritems():
            if not Prototype.prototypes[kind].has_key(key):
                Prototype.prototypes[kind][key] = value
            else:
                if Prototype.prototypes[kind][key] != value:
                    # print "old: value=" + str(Prototype.prototypes[kind][key]) + " type=" + str(type(Prototype.prototypes[kind][key]))
                    # print "new: value=" + str(value) + " type=" + str(type(value))
                    warn ("overriding prototype property "+kind+"/"+key+" to "+str(value)+" (was: "+str(Prototype.prototypes[kind][key])+")")
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
Prototype('G', radius=.05, color='gray')

Prototype('C', radius=.17, color='teal', atomic_number=6)
Prototype('C2', radius=.17, color='teal')
Prototype('CB', radius=.17, color='teal')
Prototype('H', radius=.12, color='white', atomic_number=1)
Prototype('HB', radius=.12, color='white')
Prototype('O', radius=.152, color='red')
Prototype('O1', radius=.152, color='red')
Prototype('F', radius=.135, color='pink')
Prototype('N', radius=.155, color='blue')
Prototype('P', radius=.18, color='orange')
Prototype('SI', radius=.210, color='yellow')
Prototype('Si', radius=.210, color='yellow')
Prototype('Si1', radius=.210, color='yellow')

Prototype('Si-cluster', radius=.410, color='yellow')
Prototype('particle', radius=.27, color='teal', atomic_number=6)
