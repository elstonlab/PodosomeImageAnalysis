import pickle
from mayavi.mlab import *
from mayavi import mlab
# import moviepy.editor as mpy
import os
import numpy as np



out_path = 'MayaviFigure'
out_path = os.path.abspath(out_path)
prefix = 'ani'
ext = '.png'

# mlab.options.offscreen = True


pod = pickle.load(open('stacked_pod_object_cart_slice.p', 'rb'))

pod = pod[:,:len(pod[0])//2,::-1]

fig = figure(size=(500, 500), bgcolor=(0,0,0))

surf1=contour3d(pod.T,contours=[0.25,0.4],vmin = 0.2, vmax =0.85,transparent=True,figure=fig)
surf2=contour3d(pod.T,contours=[0.6,0.8],vmin = 0.2, vmax =0.85,transparent=False,figure=fig)

cb = mlab.colorbar(nb_labels=5,label_fmt="%.2f",orientation='horizontal')
# cb.scalar_bar_representation.position = [0.6, 0.4]
# cb.scalar_bar_representation.position = [0.4, 0.4]
# # cb.scalar_bar_representation.position2 = [0.8, 0.05]

mlab.view(azimuth=90,elevation=90,distance = 320)

oa = mlab.orientation_axes(xlabel='',ylabel='',zlabel='')
# # oa.marker.set_viewport(0.0, 0.0, 0.5, 0.5)
oa.axes.normalized_label_position = [2,5,2]

cb.scalar_bar.unconstrained_font_size = True
cb.label_text_property.font_size=100

mlab.show()