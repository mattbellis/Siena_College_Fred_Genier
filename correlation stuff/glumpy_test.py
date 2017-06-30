'''
import numpy as np
from glumpy.api.matplotlib import *


fig = Figure()

left = fig.add_axes( xscale = LinearScale(clamp=True),
                        yscale = LinearScale(clamp=True),
                        zscale = LinearScale(clamp=True),
                        interface = Trackball(name="trackball"),
                        facecolor=(1,0,0,0.25), aspect=1 )

collection = PointCollection("agg")

# Add a view of the collection on the left subplot
left.add_drawable(collection)

# Add some points
collection.append(data)

fig.show()
'''

import numpy as np
from glumpy.api.matplotlib import *

data=np.loadtxt('200k_pts_xyz_data').transpose()

# Create a new figure
figure = Figure((12,12))

# Create a subplot on left, using trackball interface (3d)
left = figure.add_axes( [0.5, 0.5, 1, 1],
                        xscale = LinearScale(domain=[-1500,1500]),
                        yscale = LinearScale(domain=[-1500,1500]),
                        zscale = LinearScale(domain=[-1500,1500]),
                        interface = Trackball(name="trackball"),
                        facecolor=(0,0,0,0), aspect=1,
                        )

# Create a new collection of points
collection = PointCollection("agg",color=(0,0,0,0.25))

# Add a view of the collection on the left subplot
left.add_drawable(collection)

# Add a view of the collection on the right subplot
#right.add_drawable(collection)

# Add some points
collection.append(data)

# Show figure
figure.show()