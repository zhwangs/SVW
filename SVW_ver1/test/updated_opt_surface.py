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

Omega = 3e14
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
N_max=2

Alpha_list=np.array([])
# n loop 
#for i in range(0,len(Nested_array)):

alpha_max=get_N_SphericalWaves_given_L_max(N_max)

Row_vec=np.zeros(  (3*totaL_len, alpha_max+2)  )*1j  # I have just made 'alpha_max' to be at the column


from progressbar import ProgressBar
pbar = ProgressBar()
for n in pbar(range(1,N_max+1)):
    if n==1:
       # pass
        Row_vec[:,0]=np.ones(totaL_len*3)#+np.ones(totaL_len*3)*1j
        #Row_vec[:,1]=np.ones(totaL_len*3)*1j
    #if n==2:
    #    Row_vec[:,1]=np.ones(totaL_len*3)*1j

    for m in pbar(range(-n,n+1)):
        #print(m)
        #print(m)
        MN=get_M_N_waves(r,theta,phi, m,n, Omega, Eps, Mu, kind=1, return_Cartesian=True)
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

#%% Linear fit 
E_x=np.reshape(E_x,(totaL_len,) ,order='F')
E_y=np.reshape(E_y,(totaL_len,),order='F' )
E_z=np.reshape(E_z,(totaL_len,) ,order='F')

# point matrix 
E_plane_wave=np.vstack([E_x,E_y,E_z])

E_plane_wave=(np.reshape(E_plane_wave, (len(E_x)+len(E_y)+len(E_z),1 ), order='F'))
# ar;;;;;;;;;;;;;;;;;;;;;;;;;;;;ray of c 
E_plane_wave_ref=np.vstack([E_x,E_y,E_z])
E_plane_wave_ref_ex=(np.reshape(E_plane_wave_ref, (len(E_x)+len(E_y)+len(E_z),1 ), order='F'))
E_p_real=E_plane_wave_ref_ex.real
E_p_imag=E_plane_wave_ref_ex.imag

#%%
E_basis_real=E_plane_wave.real
# array of d
E_basis_imag=E_plane_wave.imag

E_basis_total_real_part=np.vstack((E_basis_real,E_basis_imag))


E_basis_total_imag_part=np.vstack((E_basis_imag,E_basis_real))




# Real fit 
Row_vec_real_real_part=Row_vec.real
Row_vec_imag_real_part=Row_vec.imag

Total_real_part=np.vstack((Row_vec_real_real_part,Row_vec_imag_real_part))

Fit = np.linalg.lstsq(Total_real_part, E_basis_total_real_part) 
Err=Fit[1]
C_real=Fit[0]

Fit_imag= np.linalg.lstsq(Total_real_part, E_basis_total_imag_part) 
Err=Fit_imag[1]
C_imag=Fit_imag[0]


E_cal_Creal=np.matmul(Total_real_part, C_real)#np.array([S[0]]).T) # I removed the transpose because I changed the 'alpha' index to be a column index
leng_cal=len(E_cal_Creal)
E_cal_imag_real=E_cal_Creal[int(leng_cal/2):]
E_cal_real_real=E_cal_Creal[0:int(leng_cal/2)] 
E_calculate_real=E_cal_real_real+E_cal_imag_real*1j

E_cal_Cimag=np.matmul(Total_real_part, C_imag)#np.array([S[0]]).T) # I removed the transpose because I changed the 'alpha' index to be a column index
leng_cal=len(E_cal_Cimag)
E_cal_imag_imag=E_cal_Cimag[int(leng_cal/2):]
E_cal_real_imag=E_cal_Cimag[0:int(leng_cal/2)] 
E_calculate_imag=E_cal_real_imag+E_cal_imag_imag*1j

#plane wave rec


E_p_rec=E_cal_real_real+E_cal_real_imag*1j


#E_cal_rec=np.reshape(E_calculate, (int(len(E_calculate)/3),3 ), order='F')
#E_plane_wave_rec=np.reshape(E_p_rec, (int(len(E_p_rec)/3),3 ), order='F')


#%%Test
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
f.suptitle('Plane-Wave reco (geo-size=10e-7),k-dir: [0,0,1],L_max=50', fontsize=30)

ax_1.plot(place_hod,E_cal_real_real)
ax_1.set_title('Calculated Real')
ax_2.plot(place_hod,E_p_real)
ax_2.set_title('Plane-Wave Real')
ax_6.plot(place_hod,E_cal_real_real-E_p_real)
ax_6.set_title('Difference,Real')

ax_3.plot(place_hod,(E_cal_real_imag))
ax_3.set_title('Calculated imag')
ax_4.plot(place_hod,E_p_imag)
ax_4.set_title('Plane-Wave imag')
ax_5.plot(place_hod,(E_cal_real_imag-E_p_imag))
ax_5.set_title('Difference, imag')


 #%% Linear 
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
f.suptitle('Plane-Wave reco (geo-size=10e-7),k-dir: [0,0,1],L_max=50', fontsize=30)

ax_1.plot(place_hod,E_cal_real_real)
ax_1.set_title('Calculated Real')
ax_2.plot(place_hod,E_p_real)
ax_2.set_title('Plane-Wave Real')
ax_6.plot(place_hod,E_cal_real_real-E_p_real)
ax_6.set_title('Difference,Real')

ax_3.plot(place_hod,(E_cal_real_imag))
ax_3.set_title('Calculated imag')
ax_4.plot(place_hod,E_p_imag)
ax_4.set_title('Plane-Wave imag')
ax_5.plot(place_hod,(E_cal_real_imag-E_p_imag))
ax_5.set_title('Difference, imag')


# %% Recons

C_mat=np.hstack([C_real,C_imag])
#C=C_real-C_imag
K=np.matmul(Row_vec, C_mat)
#E_cal_C=np.matmul(Row_vec, C)
Fit_C_mat= np.linalg.lstsq(K, E_p_rec) 
A=Fit_C_mat[0]

C=A[0]*C_real+A[1]*C_imag

E_cal_C=np.matmul(Row_vec, C)
from scipy.optimize import minimize

def cost_fuc(a):
    C=a[0]*C_real+a[1]*C_imag
    E_cal_C=np.matmul(Row_vec, C)
    Cost_precent=np.abs(np.sum(E_cal_C-E_plane_wave_ref_ex))
    print(Cost_precent)
    return Cost_precent

x0=np.array([1,1])
for i in range(0,100):
    cost_fuc(np.random.rand(2,1))
res_wrapped = minimize(cost_fuc,x0,method='nelder-mead')
a=res_wrapped.x
C=a[0]*C_real+a[1]*C_imag
E_cal_C=np.matmul(Row_vec, C)

#E_cal_complex=np.matmul(Total_real_part, res.x)
#%%
plt.quiver(x,y,z,E_cal_rec[:,0] ,E_cal_rec[:,1],E_cal_rec[:,2])
plt.plot(E_cal,E_basis_total_real_part)
#%% Pytorch
plt.quiver(x,y,z,E_plane_wave_rec[:,0] ,E_plane_wave_rec[:,1],E_plane_wave_rec[:,2])

#%%
f,axs=plt.subplots(2,2)
place_hod=np.linspace(0,1,int(leng_cal))
#plt.plot(place_hod,E_cal)
#plt.plot(place_hod,E_basis_total_real_part)
ax_1=axs[0,0]
ax_2=axs[0,1]

ax_1.plot(place_hod,E_cal)
ax_2.plot(place_hod,E_basis_total_real_part)

#%%

'''
import torch
import torch.nn.functional as F
# array of c 
E_basis_real=E_basis.real
# array of d
E_basis_imag=E_basis.imag

E_basis_total_real_part= np.vstack((E_basis_real,E_basis_imag))


Row_vec_real_real_part=Row_vec.real
Row_vec_imag_real_part=Row_vec.imag

Total_real_part=np.vstack((Row_vec_real_real_part,Row_vec_imag_real_part))



y=torch.from_numpy(E_basis_total_real_part)


class Net(torch.nn.Module):
    def __init__(self,n_feature,n_hidden,n_output):
        super(Net,self).__init__()

        self.hidden=torch.nn.Linear(n_feature,n_hidden)
        self.predict=torch.nn.Linear(n_hidden,n_output)

    def forward(self,x):
        x=F.relu(self.hidden(x))
        x=self.predict(x)
        return x 

net=Net(126,1,1)
print(net)

def cost_fuc(x_predicted,y):
    leng_x=len(x_predicted)
    x_predicted=x_predicted.numpy()
    y=y.numpy()
    y_predicted=np.matmul(Total_real_part,x_predicted).T 
    print(y_predicted.shape)
    print(y[:,0].shape)
    Cost_precent=2*(y_predicted-y[:,0])/(y_predicted+y[:,0])
    
    return torch.from_numpy(np.absolute(Cost_precent))


x_p=torch.linspace(0,1,Total_real_part.shape[1])
optimizer=torch.optim.SGD(net.parameters(),lr=0.2)
ls=cost_fuc(x_p,y)
pred_fuc=net(x_p)
optimizer.zero_grad()
ls.backward()
#optimizer.step()

'''












# %% using np.linalg.lstsq ------------------------------------------------
#Cost function 


#_test=get_M_N_waves(r,theta,phi, 0,1, Omega, Eps, Mu, kind=1, return_Cartesian=True)[0]
#E_basis=np.reshape(E_test, (len(x)+len(y)+len(z),1), order='F')

# array of c 
E_basis_real=E_basis.real
# array of d
E_basis_imag=E_basis.imag

E_basis_total_real_part=np.vstack((E_basis_real,E_basis_imag))

# Real fit 
Row_vec_real_real_part=Row_vec.real
Row_vec_imag_real_part=Row_vec.imag

Total_real_part=np.vstack((Row_vec_real_real_part,Row_vec_imag_real_part))


#one_array=np.ones((len(E_basis_total_real_part),1))
#Total_real_part=np.hstack([Total_real_part,one_array])



def cost_fuc(x):
    leng_x=len(x)
    #new_x=np.reshape(x,(int(leng_x/2),2))
   # real_x=x[0:int(leng_x)]
   # img_x=new_x[:,1]*1j
   # reco_x=new_x[:,0]+img_x
    E_cal_complex=np.matmul(Total_real_part, x) 
     
    Cost_precent=2*(E_cal_complex-E_basis_total_real_part.T[0,:])/(E_cal_complex+E_basis_total_real_part.T[0,:])

    return np.sum(np.absolute(Cost_precent))

initial_guess=np.ones(Row_vec.shape[1])*10**(-5)
#initial_guess=np.hstack([1,initial_guess])

E_cal_complex=np.matmul(Total_real_part, initial_guess) 
Cost_precent=2*(E_cal_complex-E_basis_total_real_part.T[0,:])/(E_cal_complex+E_basis_total_real_part.T[0,:])

res = minimize(cost_fuc, initial_guess, method='nelder-mead',
               options={'xatol': 1e-8, 'disp': True})




k=0
while k<10:
    
    c=cost_fuc(initial_guess)
    initial_guess=np.random.rand(Row_vec.shape[1])
    k+=1
    print(c)



#+np.random.rand(Row_vec.shape[1])*1j

















# %% using np.linalg.lstsq





# electric field 


#para=60
#Total_real_part=Total_real_part[np.ix_(range(0,para),range(0,len(E_basis_total_real_part)))]

#Total_real_part=Total_real_part[np.ix_(range(0,para),range(0,para))]

#Comp=np.linalg.solve(Total_real_part,E_basis_total_real_part)
#E=np.matmul(Total_real_part, Comp) 

Fit = np.linalg.lstsq(Total_real_part, E_basis_total_real_part) 
Err=Fit[1]
C=Fit[0]

#E_cal=np.matmul(Row_vec.transpose(), C)
E_cal=np.matmul(Total_real_part, C) # I removed the transpose because I changed the 'alpha' index to be a column index

# Recombied 
E_cal_complex=np.reshape(E_cal,(int(len(E_cal)/2),2))


# Precent Difference 
Err=2*(E_cal_complex-E_basis)/(E_cal_complex+E_basis)

E_summ=np.hstack((E_cal,E_basis_total_real_part))

print(E_summ)

# %% using scipy.optimize.least_squares 
#def get_fields_from_C(x, Row_vec):
import pandas as pd
pand=pd. DataFrame(E_summ)


# %%

def cost_fuc(z):
    E_cal_complex=np.matmul(Total_real_part, z) 
    Cost_precent=2*(E_cal_complex-E_basis_total_real_part.T[0,:])/(E_cal_complex+E_basis_total_real_part.T[0,:])

    return np.sum(np.absolute(Cost_precent))
def f_wrap(x):
    fx = cost_fuc(x[0] + 1j*x[1])
    return np.array([fx.real, fx.imag])

initial_guess_x=np.random.rand(Row_vec.shape[1])
initial_guess_y=np.random.rand(Row_vec.shape[1])
res_wrapped = least_squares(f_wrap, (initial_guess_x,initial_guess_y))

z = res_wrapped.x[0] + res_wrapped.x[1]*1j

# %%

 cost_fuc(initial_guess_x + 1j*initial_guess_y)

# %%
