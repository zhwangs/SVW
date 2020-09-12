
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


#%% Check M and N (Checked)
m=-1
n=2
kk=1
K=kk*np.array([1,0,1])
radius_tensor=torch.tensor(1.1,requires_grad=True)
theta_tensor=torch.tensor(0.47,requires_grad=True)
phi_tensor=torch.tensor(1.2,requires_grad=True)

S=MN_mn(radius_tensor,theta_tensor,phi_tensor,K,m,n,RG=True)
M_mn=np.conjugate(S[0]*(-1)**m)
N_mn=np.conjugate(S[1]*(-1)**m)

S2=MN_mn(radius_tensor,theta_tensor,phi_tensor,K,-m,n,RG=True)
M_mn2=S2[0] 
N_mn2=S2[1] 

print(M_mn2)
print(M_mn)
print('----')

print(N_mn2)
print(N_mn)
#%%
