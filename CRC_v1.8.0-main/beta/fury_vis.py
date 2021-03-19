
from pyMCDS_cells import pyMCDS_cells
import numpy as np
from fury import window, actor, ui

#mcds = pyMCDS_cells('output00000001.xml','data')
#mcds = pyMCDS_cells('output00000001.xml','.')
#mcds = pyMCDS_cells('output00000010.xml','.')
#mcds = pyMCDS_cells('output00000000.xml','.')
mcds = pyMCDS_cells('output00000001.xml','.')

#mcds = pyMCDS_cells('output00000001.xml','.') #  23123 cells
#mcds = pyMCDS_cells('output00000246.xml','.')  # 116038 cells
tmins = mcds.get_time()
print('time (mins)=',tmins)
print('time (days)=',tmins/1440.)

print(mcds.data['discrete_cells'].keys())
#Out[7]: dict_keys(['ID', 'position_x', 'position_y', 'position_z', 'total_volume', 'cell_type', 'cycle_model', 'current_phase', 'elapsed_time_in_phase', 'nuclear_volume', 'cytoplasmic_volume', 'fluid_fraction', 'calcified_fraction', 'orientation_x', 'orientation_y', 'orientation_z', 'polarity', 'migration_speed', 'motility_vector_x', 'motility_vector_y', 'motility_vector_z', 'migration_bias', 'motility_bias_direction_x', 'motility_bias_direction_y', 'motility_bias_direction_z', 'persistence_time', 'motility_reserved', 'oncoprotein', 'elastic_coefficient', 'kill_rate', 'attachment_lifetime', 'attachment_rate'])

# http://www.mathcancer.org/blog/paraview-for-physicell-part-1/

# The following lines assign an integer to represent
# a color, defined in a Color Map.
# sval = 0   # immune cells are yellow?
# if val[5,idx] == 1:  # [5]=cell_type
#   sval = 1   # lime green
# if (val[6,idx] == 6) or (val[6,idx] == 7):
#   sval = 0
# if val[7,idx] == 100:  # [7]=current_phase
#   sval = 3   # apoptotic: red
# if val[7,idx] > 100 and val[7,idx] < 104:
#   sval = 2   # necrotic: brownish

ncells = len(mcds.data['discrete_cells']['ID'])
print('num cells originally = ',ncells)

#xyz = np.empty((ncells,3))
xyz = np.zeros((ncells,3))
xvals = mcds.data['discrete_cells']['position_x']
yvals = mcds.data['discrete_cells']['position_y']
zvals = mcds.data['discrete_cells']['position_z']

# lets just extract half of the spheroid of tumor cells
# idx_keep = np.where(zvals < 0.0)
# ncells = len(idx_keep[0])
# print("num cells after crop = ", ncells)
# xvals = xvals[idx_keep]
# yvals = yvals[idx_keep]
# zvals = zvals[idx_keep]

xyz =np.transpose(np.array([xvals,yvals,zvals]))

# sphere V = 4/3 * pi * r^3
# r3 = V * 0.75 / pi
# r = np.cbrt(r3)
cell_radii = mcds.data['discrete_cells']['total_volume'] * 0.75 / np.pi
cell_radii = np.cbrt(cell_radii)
# cell_radii = cell_radii[idx_keep]

cell_type = mcds.data['discrete_cells']['cell_type']
# cell_type = cell_type[idx_keep]
print('cell_type min, max= ',cell_type.min(),cell_type.max())
#print(cell_type)
#cd8 = np.where(cell_type == 3.0)
#print('# cd8, macrophage, neutrophil = ',len(cd8[0]), len(macrophage[0]), len(neutrophil[0]) )

# Loop over all output files and store times and counts of cell types
#num_cd8 = np.zeros(n)
#num_neut = np.zeros(n)

rgb = np.zeros((ncells,3))
# default color: yellow
rgb[:,0] = 1
rgb[:,1] = 1
rgb[:,2] = 0
cell_phase = mcds.data['discrete_cells']['current_phase']
# cell_phase = cell_phase[idx_keep]

cycle_model = mcds.data['discrete_cells']['cycle_model']
# cycle_model = cycle_model[idx_keep]

cell_type = mcds.data['discrete_cells']['cell_type']
# cell_type = cell_type[idx_keep]

# onco = mcds.data['discrete_cells']['oncoprotein']
# onco = onco[idx_keep]
# onco_min = onco.min()
# print('onco min, max= ',onco.min(),onco.max())
# onco_range = onco.max() - onco.min()

print('cell_phase min, max= ',cell_phase.min(),cell_phase.max())  # e.g., 14.0 100.0

# This coloring is only approximately correct, but at least it shows variation in cell colors
for idx in range(ncells):
    if cell_type[idx] == 1:  # fibroblasts
        rgb[idx,0] = 0
        rgb[idx,1] = 1
        rgb[idx,2] = 1
    elif cell_type[idx] == 2:
        rgb[idx,0] = 0
        rgb[idx,1] = 1
        rgb[idx,2] = 0
    elif cell_type[idx] == 3:
        rgb[idx,0] = 1
        rgb[idx,1] = 1
        rgb[idx,2] = 0
    elif cell_type[idx] == 4:
        rgb[idx,0] = 1
        rgb[idx,1] = 0
        rgb[idx,2] = 0
        # self.yval1 = np.array( [(np.count_nonzero((mcds[idx].data['discrete_cells']['cell_type'] == 1) & (mcds[idx].data['discrete_cells']['cycle_model'] < 100) == True)) for idx in range(ds_count)] )
    # if cycle_model[idx] < 100:
    #     # rgb[idx,0] = 0.5
    #     # rgb[idx,1] = 0.5
    #     # rgb[idx,0] = 1.0 - (onco[idx] - onco_min)/onco_range
    #     rgb[idx,0] = 1.0
    #     # rgb[idx,1] = (onco[idx] - onco_min)/onco_range
    #     # rgb[idx,1] = rgb[idx,0]
    #     rgb[idx,1] = 0
    #     rgb[idx,2] = 0
    # elif cycle_model[idx] == 100:
    #     rgb[idx,0] = 1
    #     rgb[idx,1] = 0
    #     rgb[idx,2] = 0
    # elif cycle_model[idx] > 100:
    #     rgb[idx,0] = 0.54   # 139./255
    #     rgb[idx,1] = 0.27   # 69./255
    #     rgb[idx,2] = 0.075  # 19./255

#-----------------------------
scene = window.Scene()

# https://fury.gl/latest/reference/fury.actor.html?highlight=sphere#fury.actor.sphere
colors = (1,0,0)
#sphere_actor = actor.sphere(centers=xyz, colors=colors, radii=1.0)
#sphere_actor = actor.sphere(centers=xyz, colors=colors, radii=cell_radii)
sphere_actor = actor.sphere(centers=xyz, colors=rgb, radii=cell_radii)
scene.add(sphere_actor)
showm = window.ShowManager(scene,
                           size=(800, 800), reset_camera=True,
                           order_transparent=False)
showm.initialize()
showm.start()

## window.record(showm.scene, size=(900, 768), out_path="viz_timer.png")
