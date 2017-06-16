import numpy as np

from vispy import scene, visuals
from vispy.color import Color


data=np.loadtxt('200k_pts_xyz_data')


canvas = scene.SceneCanvas(keys='interactive', size=(800, 600), show=True)

# Set up a viewbox to display the cube with interactive arcball
view = canvas.central_widget.add_view()
view.bgcolor = '#000000'
view.camera = 'arcball'
view.padding = 100


markers = scene.visuals.Markers(pos=data.transpose(), face_color = "white", edge_color="white", 
                                symbol = 'disc', parent = view.scene,scaling=True)


canvas.app.run()