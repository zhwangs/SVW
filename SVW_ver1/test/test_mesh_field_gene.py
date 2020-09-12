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
gmsh_file = path_this_script_splitted[0]+'/cache/sample_mesh/mesh_4.msh'

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
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.scatter3D(Data_x, Data_y, Data_z, color='black',alpha=0.01)

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

#f = plt.figure(figsize=(15,30))
#ax = f.add_subplot(211)
#ax2 = f.add_subplot(212)

#ax.scatter(s_theta,k_theta[0:len(k_theta)-1])

#ax2.scatter(s_phi,k_phi[0:len(k_phi)-1])


# %%
# renormalize to NxN matrix 
K=s.L_grid(20)
K_theta=K[3]
K_phi=K[4]
K_radius=K[2]
R=K[1]
K_dir=K[0]
Rec=s.generate_rec_mesh(K_theta,K_phi,K_radius)
s.mesh_save_rec(path_this_script_splitted[0]+'/test/ShoMesh_data/shpmesh.mat')


# %%
Eps_r=1; Mu_r=1

Mu_vacuum = 1.25663706212e-6    # vacuum magnetic permeability
Eps_vacuum = 8.8541878128e-12   # vacuum electric permittivity
ZVAC = 376.73031346177          # impedance of free space
c = 299792458                   # the speed of light, m/s
pi = np.pi

Omega = 3e14
Eps=Eps_vacuum
Mu = Mu_vacuum
lMax = 3
polar_psi=0
polar_chi=0
a=0
b=0
E0=1
k=Omega/c#2*pi/span 



# %%
Eps_r=1; Mu_r=1

Mu_vacuum = 1.25663706212e-6    # vacuum magnetic permeability
Eps_vacuum = 8.8541878128e-12   # vacuum electric permittivity
ZVAC = 376.73031346177          # impedance of free space
c = 299792458                   # the speed of light, m/s
pi = np.pi

Omega = 3e14
Eps=Eps_vacuum
Mu = Mu_vacuum
lMax = 3
polar_psi=0#np.pi/2#np.pi/3
polar_chi=0#2*np.pi/3
a=0#np.pi/3
b=0#np.pi/3

E0=1
k=Omega/c#2*pi/span 

Geo_size=3e-4
c = 299792458 
Geo_equ_omega=2*np.pi*c/Geo_size


# %%
e=1
mu=1

#x=x.reshape((init_len**2,))
#y=y.reshape((init_len**2,))
#z=x.reshape((init_len**2,))
#print(x.shape)
#loc_array=np.array([x,y,z]).transpose()
#print(loc_array.shape)
#x=np.linspace(0,20,20)
#y=np.linspace(0,20,20)
#z=np.linspace(0,20,20)
#x, y = np.meshgrid(x, y, indexing='ij')
def plane_wave(x,y,z,polar_psi,polar_chi,a,b,E0,k):
    init_len=len(x)

    polar_chi=polar_chi*np.ones((init_len,init_len))
    polar_psi=polar_psi*np.ones((init_len,init_len))
    a=a*np.ones((init_len,init_len))
    b=b*np.ones((init_len,init_len))
    E0=E0*np.ones((init_len,init_len))
    k=k*np.ones((init_len,init_len)) 


    E_a=E0*(np.cos(polar_chi)*np.sin(polar_psi)-1j*np.sin(polar_chi)*np.cos(polar_psi))
    E_b=E0*(np.cos(polar_chi)*np.cos(polar_psi)+1j*np.sin(polar_chi)*np.sin(polar_psi))
    e_a=np.array([-np.sin(a),np.cos(a),0])
    e_b=np.array([np.cos(a)*np.cos(b),np.sin(a)*np.cos(b),-np.sin(b)])
    e_k=np.array([np.cos(a)*np.sin(b),np.sin(a)*np.sin(b),np.cos(b)])

    Dot_product=e_k[0,:,:]*x+e_k[1,:,:]*y+e_k[2,:,:]*z
    E_x=(E_a*e_a[0]+E_b*e_b[0])*np.exp(1j*Dot_product)
    E_y=(E_a*e_a[1]+E_b*e_b[1])*np.exp(1j*Dot_product)
    E_z=(E_a*e_a[2]+E_b*e_b[2])*np.exp(1j*Dot_product)



    return E_x,E_y,E_z


#%%
x,y,z=Geo_size*Rec[0],Geo_size*Rec[1],Geo_size*Rec[2]
S=plane_wave(x,y,z,polar_psi,polar_chi,a,b,E0,k)
E_x=S[0]
E_y=S[1]
E_z=S[2]

