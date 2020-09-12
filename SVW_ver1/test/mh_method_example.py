# %% Importing
#########################################################
#########################################################

import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
# any location 
Class_full_path='/home/zhw/SphEM-1/test/class_path'
sys.path.append(Class_full_path)
from mh_method_class import * 

#%% def the paths and the class 
root_path='/home/melzouka/elzouka_codes_library/Python/SphEM/'
csv_path='/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/'
processed_data_dir='/home/zhw/SphEM-1/Processed_data'
R=pointer(root_path,csv_path,processed_data_dir)
# %% define the para that can play with 
B=R.meshsimu()
geom_id=['wireSphEnds/Lx_0.3/D_0.04/']#,'wire/Lx_0.4/D_0.08/','spheroid/rx_0.04/rz_0.08/']
mesh_id=[]#['wireSphEnds/Lx_0.3/D_0.04//clscale_1_clmax_1','wireSphEnds/Lx_0.3/D_0.04//clscale_0.95_clmax_1.2','spheroid2/rx_0.04/rz_0.08//clscale_0.95_clmax_1','adhihadih']
simu_id=[]
material_id=['CONST_EPS_1.8//Au']
# S gives the pandas array with full path location and desired 4-id match for all simulation type.
S=R.file_look(geom_id,mesh_id,material_id,simu_id)
# AA gives all properties for all simulation type in three pandas arrays. 
AA=R.full_properties_look()

# choose a subset of interested properties. 
neq_prop=[]
scatter_prop=['scatter_CrossSection_m2','freq_radpers']
tmatrix_prop=[]
properties_total=[neq_prop,scatter_prop,tmatrix_prop]
# SS gives the desired properties (subset) from AA. 
SS=R.properties_look(properties_total )
# one can use df.stack() to group simulation id for better vis. 

# this gives the pandas array 
#['simulation_type', 'mesh_id', 'material_surr', 'material_particle', 'tolerance', 'RelOrient_id']
KK=R.group_by('mesh_id')
PP=R.get_prop_group_by(prop=properties_total)
SSS=R.data_see('scatterCONST_EPS_1.8Au1e-05To be codedTo be codednan',
                'scatter',scatter_prop)
# End 




# %%
