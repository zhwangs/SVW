
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



#%% Generate MN function 
m=1
n=2
kk=0.01
K=kk*np.array([0,1,1]) 
nu=torch.tensor([0.0],requires_grad=True)
x=nu2x(nu)

W_1=WignerD_fuc(x,n,0,m)

 
 

# %%
