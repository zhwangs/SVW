
# %% Importing

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


import pandas as pd
import numpy as np
all_data_orig=pd.read_pickle('/home/melzouka/elzouka_codes_library/Python/SphEM/results/digested_for_Steve/all_data_tmatrix_for_Steve_fixed_tmatrix.pkl')
Col=all_data_orig.columns

Mu_vacuum = 1.25663706212e-6    # vacuum magnetic permeability
Eps_vacuum = 8.8541878128e-12   # vacuum electric permittivity

Dir_const=1.8

c = 1/np.sqrt(Mu_vacuum*Eps_vacuum*Dir_const)
# Check L max
#L_max=np.array([all_data_orig['Lmax']])

# Unique 
freq=np.array([all_data_orig['freq_radpers']])
unique_freq=np.unique(freq,return_index=True)
leng_uniq_freq=len(unique_freq[0])


# Grouping by frequency 
Grouped_freq=all_data_orig.groupby('freq_radpers')


# Sort the frequency 
unique_freq_array_sort=np.sort(unique_freq[0])
arg_sort_freq=np.argsort(unique_freq[0])

Summary_cross_sec=np.array([0,0,0,0])

for r in range(0,leng_uniq_freq):
    Current_freq=Grouped_freq.get_group(unique_freq_array_sort[r])
    Current_L_max=np.array(Current_freq['Lmax'])
    Current_L_max_len=len(Current_L_max)
    L_len=len(Current_L_max)
    freq=unique_freq_array_sort[r]
    print(Current_L_max_len)
    k=freq/c

    for s in range(0,L_len):
        Tmatrix=np.array(Current_freq['Tmatrix'])[s][0]
        Tmatrix_sq=Tmatrix*np.conj(Tmatrix)
        C=np.sum(Tmatrix_sq)
        C=np.array([2*np.pi/(k**2)*(C)])
        C_e=np.trace(Tmatrix)
        C_e=np.array([-2*np.pi/(k**2)*(C_e)])
        C_obs=C_e.real-C
        


        Summary=np.array([k,int(Current_L_max[s]),C[0].real,C_obs[0].real])
        Summary_cross_sec=np.vstack([Summary_cross_sec,Summary])











#%% Plane Wave
def intg_fuc(x,y_1,y_2,cos_theta):
    return np.trapz(y_1,x=x)/np.pi,np.trapz(y_2,x=x)/np.pi # pi is used for cos=1. If cos-cos, replace pi by 2 

all_data_orig_plane_wave = pd.read_pickle('/home/melzouka/elzouka_codes_library/Python/SphEM/results/digested_for_Steve/all_data_scatter_for_Steve.pkl')


freq_p=np.array([all_data_orig_plane_wave['freq_radpers']])
unique_freq_p=np.unique(freq_p,return_counts=True)
leng_uniq_freq_p=len(unique_freq_p[0])


# Grouping by frequency 
Grouped_freq_p=all_data_orig_plane_wave.groupby('freq_radpers')


# Sort the frequency 
unique_freq_array_sort_p=np.sort(unique_freq_p[0])
#arg_sort_freq=np.argsort(unique_freq_p[0])

Summary_cross_sec_p=np.array([0,0,0,0])



for r in range(0,leng_uniq_freq_p):
    Current_freq=Grouped_freq_p.get_group(unique_freq_array_sort_p[r])

    freq=unique_freq_array_sort_p[r]

    k=freq/c
    plane_C=np.array(Current_freq['scatter_CrossSection_m2'])
    plane_C_obs=np.array(Current_freq['absorb_CrossSection_m2'])
    # max 

    
    max_C=np.max(plane_C)
    max_C_obs=np.max(plane_C_obs)

    Max_loc_check_C=Current_freq['scatter_CrossSection_m2']==max_C
    Max_loc_check_C_obs=Current_freq['absorb_CrossSection_m2']==max_C_obs

    #print(Max_loc_check_C)
    

    Max_array=Current_freq[Max_loc_check_C]
     


    dir_theta=np.array(Current_freq['dir_theta'])
    polar_theta=np.array(Current_freq['pol_theta'])

    #ws_sort_theta=dir_theta.reshape((int(len(dir_theta)/2),2))
    #ws_sort_plane_C=plane_C.reshape((int(len(dir_theta)/2),2))
    #ws_sort_plane_C_obs=plane_C_obs.reshape((int(len(dir_theta)/2),2))
    #ws_sort_polar_theta=polar_theta.reshape((int(len(dir_theta)/2),2))
    




    # sort angles and Cross secs
    sort_theta=np.sort(dir_theta)
    sort_theta_arg=np.argsort(dir_theta)
    

    plane_C_sort=plane_C[sort_theta_arg]
    
    plane_C_obs_sort=plane_C_obs[sort_theta_arg]

    # average polarz
    reduced_theta=sort_theta.reshape((int(len(dir_theta)/2),2))
    
    reduced_theta_ave=np.mean(reduced_theta,axis=1)

    plane_C_sort=plane_C_sort.reshape((int(len(dir_theta)/2),2))
    
    plane_C_sort_ave=np.mean(plane_C_sort,axis=1)
     

    plane_C_obs_sort=plane_C_obs_sort.reshape((int(len(dir_theta)/2),2))
    plane_C_obs_sort_ave=np.mean(plane_C_obs_sort,axis=1)


    cos_theta=1#np.cos(reduced_theta_ave) # Average method, I choose 1 here.  
    y_1=plane_C_sort_ave* np.abs(cos_theta)
    x=reduced_theta_ave
    y_2=plane_C_obs_sort_ave*np.abs(cos_theta)
    C,C_obs=intg_fuc(x ,y_1 ,y_2,cos_theta)
     
    plt.plot(x,y_1)

    Summary=np.array([k,0,C,C_obs])
    Summary_cross_sec_p=np.vstack([Summary_cross_sec_p,Summary])



#%% Testing 
unique_label=np.sort(np.unique(Summary_cross_sec[:,1]))
leng_uniq=len(unique_label)

f, axs = plt.subplots(2,1,figsize=(15,30))
ax_1=axs[0]


ax_2=axs[1]

for l in range(1,leng_uniq):
    check_loc=Summary_cross_sec[:,1]==unique_label[l]

    K=Summary_cross_sec[:,0][check_loc]
    C=Summary_cross_sec[:,2][check_loc]
    C_abs=np.abs(Summary_cross_sec[:,3][check_loc])
    ax_1.scatter(K,C,label=('L_max'+str(unique_label[l])))
    ax_2.scatter(K,C_abs,label=('L_max'+str(unique_label[l])))

Max_lim=np.max(Summary_cross_sec[:,2])
Min_lim=np.min(Summary_cross_sec[:,2])


Max_lim_a=np.max(Summary_cross_sec[:,2])
Min_lim_a=np.min(Summary_cross_sec[:,3])
#ax_1.set_xlim([xmin,xmax])
ax_1.set_ylim([Min_lim,Max_lim])
ax_2.set_ylim([Min_lim_a,Max_lim_a])

K_p=Summary_cross_sec_p[:,0]
C_p=Summary_cross_sec_p[:,2]
C_abs_p=Summary_cross_sec_p[:,3]
ax_1.plot(K_p,C_p,label='Plane Wave solution')

ax_2.plot(K_p,C_abs_p,label='Plane Wave solution')

ax_1.set_title('Average Scatter Cross Sections')
ax_2.set_title('Average Absorption Cross Sections')
#np.trapz(np.abs(cos_theta),x=x)
ax_1.set_ylabel('Scatter Cross section')
ax_1.set_xlabel('Wave vector')
ax_2.set_ylabel('Absorption Cross section')
ax_2.set_xlabel('Wave vector')
ax_1.legend()
ax_2.legend()
ax_1.grid(True)
ax_2.grid(True)


# %%
Summary_cross_sec_p
#%%
plane_C=np.array([all_data_orig_plane_wave['scatter_CrossSection_m2']])
plane_C_obs=np.array([all_data_orig_plane_wave['absorb_CrossSection_m2']])



plane_k=2*np.pi/np.array([all_data_orig_plane_wave['wavelength_nm']])

plane_freq=np.array([all_data_orig_plane_wave['freq_radpers']])
dir_x=np.array([all_data_orig_plane_wave['dir_x']])
dir_y=np.array([all_data_orig_plane_wave['dir_y']])
dir_z=np.array([all_data_orig_plane_wave['dir_z']])

dir_theta=np.array([all_data_orig_plane_wave['dir_theta']])
dir_phi=np.array([all_data_orig_plane_wave['dir_phi']])

pol_theta=np.array([all_data_orig_plane_wave['pol_theta']])




# matching index 
plane_k_unique=np.unique(plane_k)
plane_k_match_index=plane_k==plane_k_unique[-1]

dir_theta=dir_theta[plane_k_match_index]
dir_phi=dir_phi[plane_k_match_index]

pol_theta=pol_theta[plane_k_match_index]

# Sort based on angle
sort_theta=np.sort(dir_theta)
sort_theta_arg=np.argsort(dir_theta)
# Reduced the angle 
reduced_theta=sort_theta.reshape((int(len(dir_theta)/2),2))
reduced_theta=np.mean(reduced_theta,axis=1)

# Updated array based on sort angle
plane_C_match=plane_C[plane_k_match_index]
plane_C_obs_match=plane_C_obs[plane_k_match_index]

# Re-con based on sorted angle 
plane_C_match=plane_C_match[sort_theta_arg]
plane_C_obs_match=plane_C_obs_match[sort_theta_arg]

# REduced cross sections 
reduced_plane_C_match=plane_C_match.reshape((int(len(dir_theta)/2),2))
reduced_plane_C_match=np.mean(reduced_plane_C_match,axis=1)

reduced_plane_C_obs_match=plane_C_obs_match.reshape((int(len(dir_theta)/2),2))
reduced_plane_C_obs_match=np.mean(reduced_plane_C_obs_match,axis=1)

# Intg 


#reduced_theta=np.delete(reduced_theta,reduced_theta[0])
#reduced_theta=np.delete(reduced_theta,int(reduced_theta[-1]))


delta_reduced_theta=np.roll(reduced_theta, -1)-reduced_theta
delta_average_theta=np.mean(delta_reduced_theta[0:len(delta_reduced_theta)-1])

# Scattering cross sections 
cos_theta=np.cos(reduced_theta)
y=reduced_plane_C_match*np.abs(cos_theta)
x=reduced_theta
C_plane=intg_fuc(y,x)#/(np.pi)

# Abs cross sections 
y=reduced_plane_C_obs_match*np.abs(cos_theta)
x=reduced_theta
C_obs_plane=intg_fuc(y,x)#/(np.pi)




# Theta



#,return_index=True)

#reduced_plane_C=plane_C.reshape((int(len(plane_C[0])/2),2))
#reduced_plane_C_obs=plane_C_obs.reshape((int(len(plane_C_obs[0])/2),2))



print('Scatter')
print(C_plane)
print(C.real[0])
print('-------')
#print(C_e)
print('Absorption')
print(C_obs_plane)
print(C_obs.real[0])






dir_x=dir_x[plane_k_match_index]
dir_y=dir_y[plane_k_match_index]
dir_z=dir_z[plane_k_match_index]



# %%
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(dir_x,dir_y,dir_z)


# %%


# %%
