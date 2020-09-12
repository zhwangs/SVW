
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

#%% Check B  (Checked)
m=-1
n=2
theta_tensor=torch.tensor(0.0,requires_grad=True)
phi_tensor=torch.tensor(0.001,requires_grad=True)

#pre_fact=factorial(n-m)/factorial(m+n)*(-1)**m
CBP=CBP_mn(theta_tensor,phi_tensor,m,n,CBP_theta=True)
#B_con=np.conjugate(CBP[1])
#B_eq=pre_fact*B_con

#CBP=CBP_mn(theta_tensor,phi_tensor,-m,n,CBP_theta=True)


B_test=CBP[2]
#print(B_eq)
print(B_test)

#%% Check C  (Checked)
m=2
n=2
theta_tensor=torch.tensor(0.3,requires_grad=True)
phi_tensor=torch.tensor(0.47,requires_grad=True)

pre_fact=factorial(n-m)/factorial(m+n)*(-1)**m
CBP=CBP_mn(theta_tensor,phi_tensor,m,n,CBP_theta=True)
B_con=np.conjugate(CBP[0])
B_eq=pre_fact*B_con

CBP=CBP_mn(theta_tensor,phi_tensor,-m,n,CBP_theta=True)


B_test=CBP[0]
print(B_eq)
print(B_test)

#%% Check P  (Checked)
m=-1
n=2
theta_tensor=torch.tensor(0.3,requires_grad=True)
phi_tensor=torch.tensor(0.47,requires_grad=True)

pre_fact=factorial(n-m)/factorial(m+n)*(-1)**m
CBP=CBP_mn(theta_tensor,phi_tensor,m,n,CBP_theta=True)
B_con=np.conjugate(CBP[2])
B_eq=pre_fact*B_con

CBP=CBP_mn(theta_tensor,phi_tensor,-m,n,CBP_theta=True)


B_test=CBP[2]
print(B_eq)
print(B_test)

# %%
