

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


# %% Testing Bessel
x=torch.linspace(2,30,requires_grad=True)
#b_0=bessel(x,0,kind=1,derivative=0,tensor_form=False)
b_1=bessel(x,2,kind=2,derivative=0)
j1=b_1
x_test=x.detach().numpy()

# %%

fig = plt.figure()
fig.set_figheight(15)
fig.set_figwidth(15)

ax = fig.add_subplot(111)
ax.plot(x_test,j1)


# %%
