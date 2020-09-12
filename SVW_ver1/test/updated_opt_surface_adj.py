# Initiated by Zihang Wang, 


# %% Importing
#########################################################
#########################################################

import os
import sys

# get path of this script
#path_this_script = os.getcwd()
path_this_script = os.path.realpath(__file__)

# add the ./src/ path to the search path
path_this_script_splitted = os.path.split(path_this_script)
path_this_script_splitted = os.path.split(path_this_script_splitted[0])
path_to_src = os.path.join(path_this_script_splitted[0], 'src')
sys.path.append(path_to_src)  # I could have used sys.path.append('../src/'), but it didn't work with the debugger

from SphTools import *
from scipy.optimize import minimize
from scipy.optimize import least_squares,lsq_linear

from test_mesh_field_gene import *

# %% Input
#########################################################
#########################################################

# define the maximum order for expansion
L_max = 1

# define frequency, relative permittivity and permeability
Eps_r=1; Mu_r=1

Mu_vacuum = 1.25663706212e-6    # vacuum magnetic permeability
Eps_vacuum = 8.8541878128e-12   # vacuum electric permittivity
ZVAC = 376.73031346177          # impedance of free space
c = 299792458                   # the speed of light, m/s
pi = np.pi

Omega = 1e15#Geo_equ_omega
Eps=Eps_vacuum
Mu = Mu_vacuum
lMax = 3
k=Omega/c

span = 2*1/k; n_point_per_side = span*10
x_min, x_max = span,-span
y_min, y_max = span, -span
z_min, z_max = span, -span


# derrived values
Mu =   Mu_vacuum       # magnetic permeability of surrounding media
Eps =Eps_vacuum    # electric permittivity of surrounding media
Z_imped = ZVAC*np.sqrt(Mu/Eps)


# %% Calculations
#########################################################
#########################################################

# define the cooridnates of points
#x=np.linspace(x_min,x_max ,10)
#y=np.linspace(y_min,y_max,10)
#z=np.linspace(z_min,z_max,10)

totaL_len=len(x)**2
#x,y,z=np.meshgrid(x,y,z)

r,theta,phi = convert_coordinates_Cart2Sph(x, y, z)

# have to round since MN cannot take pi for theta 
#r=np.round(r,5)
theta=np.round(theta,5)
#phi=np.round(phi,5)

# %% Testing by using pure M wave 
MN_1=get_M_N_waves(r,theta,phi, -1,1, Omega, Eps, Mu, kind=1, return_Cartesian=True)
MN_2=get_M_N_waves(r,theta,phi, -1,5, Omega, Eps, Mu, kind=1, return_Cartesian=True)

E_plane=np.array([E_x,E_y,E_z]).T

F=(MN_1[0]*np.conj(MN_1[1]))
SUM_MN_theta=np.trapz(np.sum(F,axis=2)*np.sin(theta),theta,axis=0)
Sum_MN=np.trapz(SUM_MN_theta,phi)



# %% creating WV matrix
N_max=3

Alpha_list=np.array([])
# n loop 
#for i in range(0,len(Nested_array)):

alpha_max=get_N_SphericalWaves_given_L_max(N_max)

Row_vec=np.zeros(  (3*totaL_len, alpha_max+2)  )*1j  # I have just made 'alpha_max' to be at the column


#from progressbar import ProgressBar
#pbar = ProgressBar()
for n in (range(1,N_max+1)):
    #if n==1:
 
       #Row_vec[:,0]=np.ones(totaL_len*3)#+np.ones(totaL_len*3)*1j
        #Row_vec[:,1]=np.ones(totaL_len*3)*1j
    #if n==2:
    #    Row_vec[:,1]=np.ones(totaL_len*3)*1j

    for m in (range(-n,n+1)):
        #print(m)
        #print(m)
        MN=get_M_N_waves(r,theta,phi, m,n, Geo_equ_omega*0.02, Eps, Mu, kind=1, return_Cartesian=True)
        ############ I think this should be two outputs here
        F=(MN_1[0]*np.conj(MN_1[1]))
        SUM_MN_theta=np.trapz(np.sum(F,axis=2)*np.sin(theta),theta,axis=0)
        Sum_MN=np.trapz(SUM_MN_theta,phi)
        for p in [0,1]:
            #print(p)
            alpha = get_alpha_for_LMP(n, m, p)+1
            Alpha_list=np.append(Alpha_list,alpha)
            MN_c=MN[p]

            Col_1=MN_c[:,:,0]
            Col_2=MN_c[:,:,1]
            Col_3=MN_c[:,:,2]

            Col_1=np.reshape(Col_1,(totaL_len,),order='F')
            Col_2=np.reshape(Col_2,(totaL_len,),order='F')
            Col_3=np.reshape(Col_3,(totaL_len,),order='F')

            # point matrix 
            E_p=np.vstack([Col_1,Col_2,Col_3])

            E_p=np.reshape(E_p, (len(Col_1)+len(Col_2)+len(Col_3), ), order='F')
# I think you intended Col_3
            #print(Vec_fuc_W)
            # Wavefunction 
            Row_vec[:,alpha]=E_p # I just replaced 'alpha' index, to make it representative for columns ---- Row_vec[alpha,:]=Vec_fuc_W
            
            

#%%

#MN=get_M_N_waves(r[-1,:],np.round(theta[-1,:],4),theta[-1,:], m,n, Omega, Eps, Mu, kind=1, return_Cartesian=True)
            
            #print(alpha)

E_x_s=E_x
E_y_s=E_y
E_z_s=E_z

#%% Linear fit 
E_x=np.reshape(E_x_s,(totaL_len,) ,order='F')
E_y=np.reshape(E_y_s,(totaL_len,),order='F' )
E_z=np.reshape(E_z_s,(totaL_len,) ,order='F')

# point matrix 
E_plane_wave=np.vstack([E_x,E_y,E_z])

E_plane_wave=(np.reshape(E_plane_wave, (len(E_x)+len(E_y)+len(E_z),1 ), order='F'))
# ar;;;;;;;;;;;;;;;;;;;;;;;;;;;;ray of c 
E_plane_wave_ref=np.vstack([E_x,E_y,E_z])

E_plane_wave_ref_ex=(np.reshape(E_plane_wave_ref, (len(E_x)+len(E_y)+len(E_z),1 ), order='F'))
E_plane_wave_ref_ex_test=(np.reshape(E_plane_wave_ref, (3,len(E_x) ), order='F'))

E_p_test_x=E_plane_wave_ref_ex_test[0,:]
E_p_test_y=E_plane_wave_ref_ex_test[1,:]
E_p_test_z=E_plane_wave_ref_ex_test[2,:]
E_x_test=np.reshape(E_p_test_x,(len(x),len(x)) ,order='F')
E_y_test=np.reshape(E_p_test_y,(len(x),len(x)),order='F' )
E_z_test=np.reshape(E_p_test_z,(len(x),len(x)) ,order='F')

E_p_real=E_plane_wave_ref_ex.real
E_p_imag=E_plane_wave_ref_ex.imag

#%%
E_basis_real=E_plane_wave.real
# array of d
E_basis_imag=E_plane_wave.imag
E_basis_flat=np.vstack((E_basis_real,E_basis_imag))

# Real fit 
Row_vec_real_part=Row_vec.real
Row_vec_imag_part=Row_vec.imag

Upper_matrix=np.hstack((Row_vec_real_part,-Row_vec_imag_part))
Lower_matrix=np.hstack((Row_vec_imag_part,Row_vec_real_part))

Complex_Coff_matrix=np.vstack([Upper_matrix,Lower_matrix])

Fit = np.linalg.lstsq(Complex_Coff_matrix, E_basis_flat) 
C=Fit[0]




E_cal_Creal=np.matmul(Complex_Coff_matrix, C)#np.array([S[0]]).T) # I removed the transpose because I changed the 'alpha' index to be a column index
leng_cal=len(E_cal_Creal)
E_cal_imag_real=E_cal_Creal[int(leng_cal/2):]
E_cal_real_real=E_cal_Creal[0:int(leng_cal/2)] 
E_calculate_real=E_cal_real_real+E_cal_imag_real*1j

C_complex_array=C.reshape((int(len(C)/2),2),order='F')
C_complex=np.array([C_complex_array[:,0]+C_complex_array[:,1]*1j]).T
E_cal_from_complex=np.matmul(Row_vec, C_complex)
E_plane_wave_cal_rec=(np.reshape(E_cal_from_complex, (3,len(E_x)), order='F'))

E_p_x_r=E_plane_wave_cal_rec[0,:]
E_p_y_r=E_plane_wave_cal_rec[1,:]
E_p_z_r=E_plane_wave_cal_rec[2,:]
E_x_r=np.reshape(E_p_x_r,(len(x),len(x)) ,order='F')
E_y_r=np.reshape(E_p_y_r,(len(x),len(x)),order='F' )
E_z_r=np.reshape(E_x_r,(len(x),len(x)) ,order='F')

E_mag_r=np.abs(E_x_r)**2+np.abs(E_y_r)**2+np.abs(E_z_r)**2

#%%Test
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
f,axs=plt.subplots(2,3)
place_hod=np.linspace(0,1,int(leng_cal/2))
#plt.plot(place_hod,E_cal)
#plt.plot(place_hod,E_basis_total_real_part)
ax_1=axs[0,0]
ax_2=axs[0,1]
ax_3=axs[1,0]
ax_4=axs[1,1]
ax_5=axs[1,2]
ax_6=axs[0,2]
f.set_figheight(15)
f.set_figwidth(15)
f.suptitle('Plane-Wave reco (geo-size=10e-7),k-dir: [0,0,1],L_max='+ str(N_max), fontsize=30)

ax_1.plot(place_hod,E_cal_real_real)
ax_1.set_title('Calculated Real')
ax_2.plot(place_hod,E_p_real)
ax_2.set_title('Plane-Wave Real')
ax_6.plot(place_hod,np.abs(E_cal_real_real-E_p_real)**2)
ax_6.set_title('Difference,Real')

ax_3.plot(place_hod,(E_cal_imag_real))
ax_3.set_title('Calculated imag')
ax_4.plot(place_hod,E_p_imag)
ax_4.set_title('Plane-Wave imag')
ax_5.plot(place_hod,(np.abs(E_cal_imag_real-E_p_imag)**2))
ax_5.set_title('Difference, imag')





# %%
fig = plt.figure()
fig.suptitle('Plane-Wave reco (geo-size=10e-7),L_max='+ str(N_max), fontsize=30)
ax_1 = fig.add_subplot(221)
ax_2 = fig.add_subplot(222)
ax_3 = fig.add_subplot(223)
ax_4 = fig.add_subplot(224, projection='3d')


fig.set_figheight(15)
fig.set_figwidth(15)

ax_1.quiver(x,y,E_x_r ,E_y_r)
ax_1.set_title('x-y plane')
ax_1.set_xlabel('x (meter)')
ax_1.set_ylabel('y (meter)')

ax_2.quiver(x,z,E_x_r ,E_z_r)
ax_2.set_title('x-z plane')
ax_2.set_xlabel('x (meter)')
ax_2.set_ylabel('z (meter)')

ax_3.quiver(y,z,E_y_r ,E_z_r)
ax_3.set_title('y-z plane')
ax_3.set_xlabel('y (meter)')
ax_3.set_ylabel('z (meter)')




norm = Normalize()
colors = norm(E_mag_r**2) # auto-adjust true radius into [0,1] for color mapping

cmap = cm.get_cmap("coolwarm")
ax_4.plot_surface(x, y, z, linewidth=0, facecolors=cmap(colors), shade=True, alpha=1)

# the surface is not mappable, we need to handle the colorbar manually
mappable = cm.ScalarMappable(cmap=cmap)
mappable.set_array(E_mag_r**2)
fig.colorbar(mappable, shrink=0.2, aspect=5,label='Field strength (E^2)')
 # auto-adjust true radius into [0,1] for color mapping
ax_4.quiver(x, y, z, E_x_r, E_y_r, E_z_r, length=0.2*x.max())
ax_4.view_init(elev=15, azim=45)
ax_4.set_title('3D view, vector field')
ax_4.set_xlabel('x (meter)')
ax_4.set_ylabel('y (meter)')
ax_4.set_zlabel('z (meter)')
# %%


fig = plt.figure()
fig.set_figheight(15)
fig.set_figwidth(15)

ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=15, azim=45)#np.pi/4)

norm = Normalize()
colors = norm(E_mag_r) # auto-adjust true radius into [0,1] for color mapping

cmap = cm.get_cmap("coolwarm")
ax.plot_surface(x, y, z, linewidth=0, facecolors=cmap(colors), shade=True, alpha=1)

# the surface is not mappable, we need to handle the colorbar manually
mappable = cm.ScalarMappable(cmap=cmap)
mappable.set_array(E_mag_r**2)
fig.colorbar(mappable, shrink=0.2, aspect=5,label='Field strength (E^2)')
 # auto-adjust true radius into [0,1] for color mapping
ax.quiver(x, y, z, E_x_r, E_y_r, E_z_r, length=0.2*x.max())

ax.set_title('3D view, vector field')
ax.set_xlabel('x (meter)')
ax.set_ylabel('y (meter)')
ax.set_zlabel('z (meter)')




# %%
