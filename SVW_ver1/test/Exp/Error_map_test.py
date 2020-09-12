
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


gmsh_file = path_this_script_splitted[0]+'/cache/sample_mesh/mesh_4.msh'

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
 



# %%

def error_map(n_max, polar_psi, polar_chi,E0,K, K_mag_range,theta, phi, radius,radius_mag_range):
    '''Def error map respect to the radius and wave vector'''
    Omega = 3e8
    Error_map_x=np.zeros((len(K_mag_range),len(radius_mag_range)))
    Error_map_y=np.zeros((len(K_mag_range),len(radius_mag_range)))
    Error_map_z=np.zeros((len(K_mag_range),len(radius_mag_range)))

    theta_tensor=torch.tensor(theta, requires_grad=True)
    phi_tensor=torch.tensor(phi, requires_grad=True)

 
    for k_mag_index in range(0,len(K_mag_range)):
        
        current_K=K*K_mag_range[k_mag_index]


        for r_mag_index in range(0,len(radius_mag_range)):

            current_radius=radius*radius_mag_range[r_mag_index]
 
            Rec=s.generate_rec_mesh(theta, phi, current_radius)
            x=Rec[0]
            y=Rec[1]
            z=Rec[2]


            
            radius_tensor=torch.tensor(current_radius, requires_grad=True)
            K2=plane_wave_dir(polar_psi,polar_chi,E0,current_K)
            RB=plane_wave_grid(x,y,z,K2,current_K,omega=Omega,t0=0)
            E_x_r=RB[0] # x-y plane 
            E_y_r=RB[1]
            E_z_r=RB[2]
            E_rec=rec_plane_wave(radius_tensor,theta_tensor,phi_tensor,current_K,K2,n_max)


            E_rec_x=E_rec[0]
            E_rec_y=E_rec[1]
            E_rec_z=E_rec[2]

            Error_x=np.sum((np.abs(E_x_r.real-E_rec_x.real))**2)/np.sum(E_x_r.real+E_rec_x.real)
            Error_y=np.sum((np.abs(E_y_r.real-E_rec_y.real))**2)/np.sum(E_y_r.real+E_rec_y.real)
            Error_z=np.sum((np.abs(E_z_r.real-E_rec_z.real))**2)/np.sum(E_y_r.real+E_rec_y.real)
            print(Error_x)
            Error_map_x[k_mag_index,r_mag_index]=Error_x
            Error_map_y[k_mag_index,r_mag_index]=Error_y
            Error_map_z[k_mag_index,r_mag_index]=Error_z

    return Error_map_x,Error_map_y,Error_map_z
 



 
#RR=angle_from_ek(k)

 
 # %%
# %%
kk=0.01
K=kk*np.array([0,1,1])
n_max=6
K_mag_range=np.linspace(0.01,1,2)#np.array([0.01,0.1,1,10])
radius_mag_range=np.linspace(100,10,2)#np.array([1000,100,10,1])

polar_psi=0#np.pi/2#np.pi/3
polar_chi=np.pi/4

lbd=2*np.pi/kk
E0=1


Er=error_map(n_max, polar_psi, polar_chi,E0,K, K_mag_range,theta, phi, radius,radius_mag_range)

#%%
import seaborn as sns


#%%
ax = sns.heatmap(Er[0], linewidth=0)
plt.show()
# %%
