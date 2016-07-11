import pygmsh as pg
import numpy as np

geom = pg.Geometry()

# Draw a cross.
poly = geom.add_polygon([
    [0.0,   0.5, 0.0],
    [-0.1,  0.1, 0.0],
    [-0.5,  0.0, 0.0],
    [-0.1, -0.1, 0.0],
    [0.0,  -0.5, 0.0],
    [0.1,  -0.1, 0.0],
    [0.5,   0.0, 0.0],
    [0.1,   0.1, 0.0]
    ],
    lcar=0.05
    )

axis = [0, 0, 1]

geom.extrude(
    'Surface{%s}' % poly,
    translation_axis=axis,
    rotation_axis=axis,
    point_on_axis=[0, 0, 0],
    angle=2.0 / 6.0 * np.pi
    )

points, cells = pg.generate_mesh(geom)