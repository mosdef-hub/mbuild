import networkx as nx
import numpy as np
import numpy.ma as ma
from mayavi import mlab


# numpy array of x,y,z positions in sorted node order

pos0 = np.array([0,0,0])
pos1 = np.array([1,0,0])
pos2 = np.array([0,1,0])
pos3 = np.array([-1,0,0])
pos4 = np.array([0,-1,0])

xyz=np.array([pos0, pos1, pos2, pos3, pos4])

# scalar color
scalars=np.array([0,1,2,3,4])

mlab.figure(1, bgcolor=(0, 0, 0))
mlab.clf()

pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
                    scalars,
                    scale_factor=.1,
                    scale_mode='scalar',
                    colormap='Blues',
                    resolution=20)



conns = [(0,1), (0,2), (0,3)]
conn_scalars = np.array([1,2,3,4])
pts.mlab_source.dataset.lines = np.array(conns)



# pts.mlab_source.dataset.lines.add_array(conn_scalars)

tube = mlab.pipeline.tube(pts, tube_radius=0.1)
# tube.filter.radius_factor = 1.
tube.filter.vary_radius = 'vary_radius_by_scalar'

mlab.pipeline.surface(tube, colormap='winter')





# import pdb
# pdb.set_trace()

mlab.show() # interactive window