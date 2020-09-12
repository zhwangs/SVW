
#%%
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


rad_test=100#2000




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
#s.angle_radius_mesh()
#Est=s.G_quad_mesh_N_esti()
#k_theta=Est[0] # theta unique 
#k_phi=Est[2] # phi unique 
#k_theta_mean=Est[1] # theta mean 
#k_phi_mean=Est[3] # phi mean 

#s_theta=np.linspace(0,1,len(k_theta)-1)
#s_phi=np.linspace(0,1,len(k_phi)-1)
 
K=s.L_grid(20)
theta=K[3]
phi=K[4]
radius=K[2]*rad_test 

Rec=s.generate_rec_mesh(K[3], K[4], K[2])


 
x=Rec[0]
y=Rec[1]
z=Rec[2]

#%% Mesh 
kk=0.0001
K_mesh=K[0]
A=Emap_K(K_mesh,n_max=6,polar_psi=0,polar_chi=0,E0=1)

# %%

K=kk*np.array([1,1,1])
n_max=6
polar_psi=0#np.pi/2#np.pi/3
polar_chi=0#np.pi/4

lbd=2*np.pi/kk
E0=1
K2=plane_wave_dir(polar_psi,polar_chi,E0,K)



#%% Plane Wave 

Mu_vacuum = 1.25663706212e-6    # vacuum magnetic permeability
Eps_vacuum = 8.8541878128e-12   # vacuum electric permittivity
ZVAC = 376.73031346177          # impedance of free space
c = 3e8#299792458                   # the speed of light, m/s
pi = np.pi

Omega = 3e8
Eps=Eps_vacuum
Mu = Mu_vacuum

x=Rec[0,:] 
y=Rec[1,:] 
z=Rec[2,:] 


RB=plane_wave_grid(x,y,z,K2,K,omega=Omega,t0=0)
E_x_r=RB[0] # x-y plane 
E_y_r=RB[1]
E_z_r=RB[2]



#%% Plane wave rec 
radius_tensor=torch.tensor(radius, requires_grad=True)
theta_tensor=torch.tensor(theta, requires_grad=True)
phi_tensor=torch.tensor(phi, requires_grad=True)

n_max=6


E_rec=rec_plane_wave(radius_tensor,theta_tensor,phi_tensor,K,K2,n_max)

E_x_r_rec=E_rec[0] # x-y plane 
E_y_r_rec=E_rec[1]
E_z_r_rec=E_rec[2]


# %%
Error_x=np.sum((E_x_r_rec.real-E_x_r.real)**2/np.abs(E_x_r_rec.real+E_x_r.real))/E_x_r_rec.reshape(-1).shape[0]
Error_y=np.sum((E_y_r_rec.real-E_y_r.real)**2/np.abs(E_x_r_rec.real+E_x_r.real))/E_y_r_rec.reshape(-1).shape[0]
Error_z=np.sum((E_z_r_rec.real-E_z_r.real)**2/np.abs(E_x_r_rec.real+E_x_r.real))/E_z_r_rec.reshape(-1).shape[0]

print(Error_x,Error_y,Error_z)
# %%
