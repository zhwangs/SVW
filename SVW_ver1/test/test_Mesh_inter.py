# %% Importing
#########################################################
#########################################################

import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

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
from Mesh_inter import *



#%%
gmsh_file = '/Users/zhwang/Desktop/SphEM-Mesh_angle/cache/sample_mesh/mesh_3.msh'

# Class 
s=Mesh_inter(gmsh_file)
# interpolation (triangle mesh)
Data_g=s.mesh_interpo(20)
Data_x=Data_g[0]
Data_y=Data_g[1]
Data_z=Data_g[2]


# %%
# Orignial mesh

Org_data=Data_g[3]
fig = plt.figure()
ax = plt.axes(projection="3d")
Center_array=np.array(list(Org_data))
ax.scatter3D(Center_array[:,0],Center_array[:,1],Center_array[:,2],alpha=1)

# %% Triangle Mesh 
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(Data_x, Data_y, Data_z, color='black',alpha=0.01)

# %% Mesh Estimate in angle 
# Create theta phi mesh 
s.angle_radius_mesh()
Est=s.G_quad_mesh_N_esti()
k_theta=Est[0] # theta unique 
k_phi=Est[2] # phi unique 
k_theta_mean=Est[1] # theta mean 
k_phi_mean=Est[3] # phi mean 

s_theta=np.linspace(0,1,len(k_theta)-1)
s_phi=np.linspace(0,1,len(k_phi)-1)

f = plt.figure(figsize=(15,30))
ax = f.add_subplot(211)
ax2 = f.add_subplot(212)

ax.scatter(s_theta,k_theta[0:len(k_theta)-1])

ax2.scatter(s_phi,k_phi[0:len(k_phi)-1])


# %%
# renormalize to NxN matrix 
K=s.L_grid(50)
K_theta=K[3]
K_phi=K[4]
K_radius=K[2]
R=K[1]
K_dir=K[0]
Rec=s.generate_rec_mesh(K_theta,K_phi,K_radius)
#s.mesh_save_rec()


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(K_dir[0], K_dir[1], K_dir[2], color='black',alpha=0.2)
# %%
# NxN shape x,y,z
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(R[:,0], R[:,1], R[:,2], color='black',alpha=0.2)
# %%
# NxN shape r,theta,phi
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(Rec[0], Rec[1], Rec[2], color='black',alpha=0.2)


# %%
