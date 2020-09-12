
#%% (Checked)
import os
import sys
# get path of this script
#path_this_script = os.getcwd()
path_this_script = os.path.realpath(__file__)

# add the ./src/ path to the search path
path_this_script_splitted = os.path.split(path_this_script)
this_script_filename = path_this_script_splitted[1]
path_this_script_splitted = os.path.split(path_this_script_splitted[0])
path_to_src = os.path.join(path_this_script_splitted[0], 'test/Exp')
sys.path.append(path_to_src)  # I could have used sys.path.append('../src/'), but it didn't work with the debugger

from WD_fuc import * 


rad_test=2000




# %%


import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.special import spherical_jn, spherical_yn

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D

# get path of this script
#path_this_script = os.getcwd()
path_this_script = os.path.realpath(__file__)

# add the ./src/ path to the search path
path_this_script_splitted = os.path.split(path_this_script)
this_script_filename = path_this_script_splitted[1]
path_this_script_splitted = os.path.split(path_this_script_splitted[0])
path_this_script_splitted = os.path.split(path_this_script_splitted[0])

path_to_src = os.path.join(path_this_script_splitted[0], 'src')
sys.path.append(path_to_src)   
path_to_cache = os.path.join(path_this_script_splitted[0], 'cache')
from Mesh_inter import *

#%%  create mesh 


gmsh_file = path_this_script_splitted[0]+'/cache/sample_mesh/mesh_5.msh'

# Class 
s=Mesh_inter(gmsh_file)
# interpolation (triangle mesh)
Data_g=s.mesh_interpo(20)
Data_x=Data_g[0]
Data_y=Data_g[1]
Data_z=Data_g[2]



# Orignial mesh

Org_data=Data_g[3]
 
# Create theta phi mesh 
s.angle_radius_mesh()
Est=s.G_quad_mesh_N_esti()
k_theta=Est[0] # theta unique 
k_phi=Est[2] # phi unique 
k_theta_mean=Est[1] # theta mean 
k_phi_mean=Est[3] # phi mean 

s_theta=np.linspace(0,1,len(k_theta)-1)
s_phi=np.linspace(0,1,len(k_phi)-1)
 
K=s.L_grid(20)
theta=K[3]
phi=K[4]
radius=K[2]*rad_test 
 
Rec=s.generate_rec_mesh(theta, phi, radius)
 
x=Rec[0]
y=Rec[1]
z=Rec[2]
 
radius_tensor=torch.tensor(radius, requires_grad=True)
theta_tensor=torch.tensor(theta, requires_grad=True)
phi_tensor=torch.tensor(phi, requires_grad=True)

#%% Generate MN function 
m=0
n=1
kk=0.01
K=kk*np.array([1,0,0])

S=MN_mn(radius_tensor,theta_tensor,phi_tensor,K,m,n,RG=True)
M_mn=S[0]
N_mn=S[1]




#%% Testing with AB coeff
polar_psi=0#np.pi/2#np.pi/3
polar_chi=0#np.pi/3

lbd=2*np.pi/kk
E0=1
K2=plane_wave_dir(polar_psi,polar_chi,E0,K)


 


M_mn_x=M_mn[0,:,:]#*AB[0]
M_mn_y=M_mn[1,:,:] 
M_mn_z=M_mn[2,:,:] 

N_mn_x=N_mn[0,:,:]#*AB[1]
N_mn_y=N_mn[1,:,:] 
N_mn_z=N_mn[2,:,:] 


MN_x=M_mn_x+N_mn_x
MN_y=M_mn_y+N_mn_y
MN_z=M_mn_z+N_mn_z


#%%

fig = plt.figure()
fig.set_figheight(15)
fig.set_figwidth(15)

ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=15, azim=-45)#np.pi/4)



E_mag_r=MN_x.real#np.sqrt(np.abs(E_x_r)**2+np.abs(E_y_r)**2+np.abs(E_z_r)**2)
norm = Normalize()
colors =norm(E_mag_r) # auto-adjust true radius into [0,1] for color mapping

cmap = cm.get_cmap("coolwarm")
ax.plot_surface(x, y, z, linewidth=0, facecolors=cmap(colors), shade=False, alpha=1)

# the surface is not mappable, we need to handle the colorbar manually
mappable = cm.ScalarMappable(cmap=cmap)
mappable.set_array(E_mag_r)
fig.colorbar(mappable, shrink=0.2, aspect=5,label='Field strength (E^2)')
 # auto-adjust true radius into [0,1] for color mapping

#ax.quiver(x_r, y_r, z_r, M_mn_x.real, M_mn_y.real, M_mn_z.real, length=0.1*x.max(),linewidths=0.6*x.max())

ax.set_title('3D view, vector field')
ax.set_xlabel('x (meter)')
ax.set_ylabel('y (meter)')
ax.set_zlabel('z (meter)')

 
# %%


# %%
