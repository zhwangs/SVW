
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
from scipy.special import spherical_jn, spherical_yn

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

#%%  plot for angular function 
phi_ek=torch.linspace(0,np.pi,100,requires_grad=True)
m=2
n=3
PI=pi_mn(phi_ek,m,n).detach().numpy()
TAU=tau_mn(phi_ek,m,n).detach().numpy()

A=PI
phi_ek=phi_ek.detach().numpy()

x=A*np.cos(phi_ek)
y=A*np.sin(phi_ek)

plt.plot(x,y)

A=TAU
x=A*np.cos(phi_ek)
y=A*np.sin(phi_ek)
plt.plot(x,y)

