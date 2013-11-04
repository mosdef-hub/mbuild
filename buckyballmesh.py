__author__ = 'sallai'

from tile import *
from mbuild.ndimmesh import *
from mbuild.twodimmesh import *
from buckyballtile import *

if __name__ == "__main__":
    m = TwoDimMesh.create(BuckyBallTile.create, 41, 41,
                          left_port_name='left_male_port',
                          right_port_name='right_female_port',
                          top_port_name='top_male_port',
                          bottom_port_name='bottom_female_port')
    # print m.atoms()
    m.plot(labels=False, verbose=False)






