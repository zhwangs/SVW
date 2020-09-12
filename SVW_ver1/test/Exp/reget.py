
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







# %%


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
path_this_script_splitted = os.path.split(path_this_script_splitted[0])

path_to_src = os.path.join(path_this_script_splitted[0], 'src')
sys.path.append(path_to_src)   
path_to_cache = os.path.join(path_this_script_splitted[0], 'cache')
from Mesh_inter import *



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
#fig = plt.figure()
#ax = plt.axes(projection="3d")
#Center_array=np.array(list(Org_data))
#ax.scatter3D(Center_array[:,0],Center_array[:,1],Center_array[:,2],alpha=1)

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


# renormalize to NxN matrix 
K=s.L_grid(20)
theta=K[3]
phi=K[4]
radius=K[2]*2000
#K_dir=K[0]
Rec=s.generate_rec_mesh(theta, phi, radius)
#s.mesh_save_rec(path_this_script_splitted[0]+'/test/ShoMesh_data/shpmesh.mat')
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.scatter3D(Rec[0], Rec[1], Rec[2], color='black',alpha=0.2)



# %%
radius_tensor=torch.tensor(radius, requires_grad=True)
theta_tensor=torch.tensor(theta, requires_grad=True)
phi_tensor=torch.tensor(phi, requires_grad=True)



# %%
WD=WignerD_fuc(theta_tensor,3,2,1)

n=1
m=0
#dd=unit2cart(s,s)

R=CBP_mn(theta_tensor,phi_tensor,m,n,CBP_theta=True)


# %%
Mu_vacuum = 1.25663706212e-6    # vacuum magnetic permeability
Eps_vacuum = 8.8541878128e-12   # vacuum electric permittivity
ZVAC = 376.73031346177          # impedance of free space
c = 3e8#299792458                   # the speed of light, m/s
pi = np.pi

Omega = 3e8
Eps=Eps_vacuum
Mu = Mu_vacuum

polar_psi=0#np.pi/2#np.pi/3
polar_chi=0#np.pi/3

E0=1
k=Omega/c#2*pi/span 
#x=np.linspace(0,10,3)
#y=np.linspace(0,2*np.pi,10)
#z=np.linspace(0,2*np.pi,10)

x=Rec[0,:] 
y=Rec[1,:] 
z=Rec[2,:] 
kk=0.1
K=kk*np.array([0,1,0])
lbd=2*np.pi/kk
K2=plane_wave_dir(polar_psi,polar_chi,E0,K)
#RR=angle_from_ek(k)

RB=plane_wave_grid(x,y,z,K2,K,omega=Omega,t0=0)
E_x_r=RB[0] # x-y plane 
E_y_r=RB[1]
E_z_r=RB[2]

AB=ab_matrix(n,m,K2,K)
 # %%
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D



#plt.quiver(x, y, E_x, E_y)


fig = plt.figure()
fig.set_figheight(15)
fig.set_figwidth(15)

ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=15, azim=45)#np.pi/4)

x_r=x
y_r=y
z_r=z



E_mag_r=E_x_r.real#np.sqrt(np.abs(E_x_r)**2+np.abs(E_y_r)**2+np.abs(E_z_r)**2)
norm = Normalize()
colors =norm(E_mag_r) # auto-adjust true radius into [0,1] for color mapping

cmap = cm.get_cmap("coolwarm")
ax.plot_surface(x, y, z, linewidth=0, facecolors=cmap(colors), shade=False, alpha=1)

# the surface is not mappable, we need to handle the colorbar manually
mappable = cm.ScalarMappable(cmap=cmap)
mappable.set_array(E_mag_r)
fig.colorbar(mappable, shrink=0.2, aspect=5,label='Field strength (E^2)')
 # auto-adjust true radius into [0,1] for color mapping

#ax.quiver(x_r, y_r, z_r, E_x_r.real, E_y_r.real, E_z_r.real, length=0.1*x.max(),linewidths=0.6*x.max())

ax.set_title('3D view, vector field')
ax.set_xlabel('x (meter)')
ax.set_ylabel('y (meter)')
ax.set_zlabel('z (meter)')

#bessel(radius_tensor,n,kind=1,derivative=1,tensor_form=False)
# %%
#S=MN_mn(radius_tensor,theta_tensor,phi_tensor,K,m,n,RG=True)
# %%

EEx=E_z_r.reshape(-1,).real
xxx=y.reshape(-1,)
plt.scatter(xxx,EEx)


# %% Sphere vector wave 
n_max=7
K=np.array([0,1,0])
E_rec=rec_plane_wave(radius_tensor,theta_tensor,phi_tensor,K,K2,n_max)
#%%
# Testing 

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D



fig = plt.figure()
fig.set_figheight(15)
fig.set_figwidth(15)

ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=15, azim=-45)#np.pi/4)

x_r=x
y_r=y
z_r=z



E_mag_r=E_rec[0].real#np.sqrt(np.abs(E_x_r)**2+np.abs(E_y_r)**2+np.abs(E_z_r)**2)
norm = Normalize()
colors =(E_mag_r) # auto-adjust true radius into [0,1] for color mapping

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


#%% Plot Sphere vector wave

