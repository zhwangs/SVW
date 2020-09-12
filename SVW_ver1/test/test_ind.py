import os
import sys

# get path of this script
#path_this_script = os.getcwd()
path_this_script = os.path.realpath(__file__)

# add the ./src/ path to the search path
path_this_script_splitted = os.path.split(path_this_script)
this_script_filename = path_this_script_splitted[1]
path_this_script_splitted = os.path.split(path_this_script_splitted[0])
path_to_src = os.path.join(path_this_script_splitted[0], 'src')
sys.path.append(path_to_src)  # I could have used sys.path.append('../src/'), but it didn't work with the debugger
path_to_cache = os.path.join(path_this_script_splitted[0], 'cache')
from SphTools import *
from GeometryTools import *








# %% test get_points_from_gmsh_file
gmsh_file = '/Users/zhwang/Desktop/SphEM/cache/sample_mesh/mesh_4.msh'
mesh_points = get_points_from_gmsh_file(gmsh_file)[0]
mesh_object=get_points_from_gmsh_file(gmsh_file)[1]



# %% test get_goemetric_center
center1 = get_goemetric_center(mesh_points)
center2 = get_goemetric_center(mesh_points[...,0], mesh_points[...,1], mesh_points[...,2])

Tir_array=mesh_object.cells_dict['triangle']
mesh_points=mesh_object.points
#mesh_points=mesh_points-center1
# angle creation 
Theta=np.linspace(0,np.pi,10)
Phi=np.linspace(0,2*np.pi,15)
# angle 
theta,phi=np.meshgrid(Theta,Phi)


s=Mesh_inter(gmsh_file)
#s.tri_fuc()
# %%
s.tri_fuc()

# %%
