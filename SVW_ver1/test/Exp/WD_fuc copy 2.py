# here

'''
You need to install the following packages:

pip install numpy
pip install scipy
pip install h5py
pip install pandas

pip install matplotlib


pip install pypandoc

'''
# %% importing modules
#################################################################
#################################################################

# general ------------------------------------------------------------
import numpy as np
import scipy as sp
import h5py
import re
import time
import datetime
import os
import pathlib
from numpy.polynomial.legendre import legval, legder
from scipy.special import lpmv, lpmn, clpmn, jv, jvp, yv, yvp, hankel1, hankel2, h1vp, h2vp, factorial, sph_harm
import sys
import pandas as pd

from pathlib import Path

import json

import runpy

from time import strftime,time
from scipy.special import spherical_jn, spherical_yn

# plotting -------------------------------------------------------------
import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

font = {    
        'size'   : 5}

matplotlib.rc('font', **font)


# global definitions and constants --------------------------------------
pi = np.pi
Mu_vacuum = 1.25663706212e-6    # vacuum magnetic permeability
Eps_vacuum = 8.8541878128e-12   # vacuum electric permittivity
ZVAC = 376.73031346177          # impedance of free space
c = 299792458                   # the speed of light, m/s
SMALL_NUMBER = 1e-10 ; # SMALL_NUMBER = np.finfo(float).eps * 1e1
Stefan_Boltzmann_const = 5.670367e-8 # Stefan-Boltzmann constant


import torch 
from scipy.special import factorial
from torch.autograd import grad

def nu2x(nu):
    '''
    equation B8 gives 
    cos(1/2 nu)=sqrt(1/2+x/2), sin(1/2 nu)=sqrt(1/2-x/2)
    input has to be a tensor object

    Checked 
    '''
    nu=nu.to(torch.float)
    x=torch.cos(nu)
    #x=2*((torch.cos(nu/2))**2)-1

    return x
def x2nu(x):
    '''
    Inverse 
    equation B8 gives 
    cos(1/2 nu)=sqrt(1/2+x/2), sin(1/2 nu)=sqrt(1/2-x/2)
    input has to be a tensor object
    '''
    #x=x.to(torch.float)
    nu=2*torch.acos(torch.sqrt(1/2+x/2))
    return nu

def A_s_mn(s,m,n):
    ''' 
    equation B10. 
    s>=0, -s<=m. n<=s 
    '''    
    A_smn=((-1)**(s-n)/2**s)*np.sqrt((factorial(s+n))/(factorial(s+m)*factorial(s-m)*factorial(s-n)))
#    if s<0 or -s>m or n>s:
 #       A_smn=0
    return A_smn

 
def nth_derivative(f, wrt, n):
    '''
    Taking n'th derivative for a flatten tensor 

    '''
    
    if n>0:
        f=f.sum()
        for i in range(n):

            grads = grad(f, wrt, create_graph=True)[0]
            f = grads.sum()
    else:
        grads=f

    return grads

def  WignerD_fuc(x,s,m,n):
    
    A_smn=A_s_mn(s,m,n)
    frist_part=(1-x)**((m-n)/2)*(1+x)**(-(m+n)/2)
    second_part=(1+x)**(s+m)*(1-x)**(s-m)
    #print(frist_part)
    # check int s-n>=0
    der_num=s-n
    
    der_second_part=nth_derivative(second_part,x,n=der_num)
 

    WingerD_out=A_smn*frist_part*der_second_part*(-1)**s
    try:
        WingerD_out_test=WingerD_out.detach().numpy()

        if np.isnan(WingerD_out_test):
            if m==n:
                WingerD_out=torch.tensor(1,requires_grad=True)
            else:
                WingerD_out=torch.tensor(0,requires_grad=True)
        else:
            pass

    except:
        pass

            

    return WingerD_out

#---------SVW --------------------------


def pi_mn(nu,m,n):
    x=nu2x(nu)
    #print(WignerD_fuc(x,n,0,m))
    pi_val=(m/torch.sin(nu))*WignerD_fuc(x,n,0,m)     
    return pi_val

def tau_mn(nu,m,n):
    '''
    d/dnu=dx/dnu 1/dx'''
    
    x=nu2x(nu)
    f=WignerD_fuc(x,n,0,m)
    S=nth_derivative(f, x, 1)
    P=nth_derivative(x, nu, 1)
    return S*P


def unit2cart(nu,phi,detach=True):
    if detach:
        nu=nu.detach().numpy()
        phi=phi.detach().numpy()
    else:
        pass
    r_hat=np.array([np.sin(nu)*np.cos(phi),np.sin(nu)*np.sin(phi),np.cos(nu)])
    nu_hat=np.array([np.cos(nu)*np.cos(phi),np.cos(nu)*np.sin(phi),-np.sin(nu)])
    phi_hat=np.array([-np.sin(phi),np.cos(phi),np.zeros(phi.shape)])
    return r_hat, nu_hat, phi_hat


def CBP_mn(nu,phi,m,n,CBP_theta=False):
    '''
    C18 C19 and C20. 
    CBP_theta gives only theta dependence when false 
    '''
    x=nu2x(nu)
 
    
    pi_mn_var=pi_mn(nu,m,n).detach().numpy()

    pi_mn_var=np.array([pi_mn_var,pi_mn_var,pi_mn_var])
    tau_mn_var=tau_mn(nu,m,n).detach().numpy()
    tau_mn_var=np.array([tau_mn_var,tau_mn_var,tau_mn_var])
    

    nu_d=nu.detach().numpy()
    phi_d=phi.detach().numpy()
    unit_vec=unit2cart(nu_d,phi_d,detach=False)
    nu_hat=unit_vec[1]
    phi_hat=unit_vec[2]
    r_hat=unit_vec[0]


    WD_mn=WignerD_fuc(x,n,0,m).detach().numpy()
 
    
    WD_mn=np.array([WD_mn,WD_mn,WD_mn])

    P_factor=((-1)**m)*np.sqrt(factorial(n+m)/factorial(n-m))*np.exp(1j*m*phi_d)
    P_factor=np.array([P_factor,P_factor,P_factor])


    B_mn=tau_mn_var*nu_hat+1j*pi_mn_var*phi_hat
    C_mn=-tau_mn_var*phi_hat+1j*pi_mn_var*nu_hat

    B_mn=B_mn
    C_mn=C_mn
    P_mn=r_hat*WD_mn
    if CBP_theta==True:

        B_mn=B_mn*P_factor
        C_mn=C_mn*P_factor

        P_mn=P_mn*P_factor
    else:
        pass 




 
    return  C_mn,B_mn,P_mn


def gamma_mn(m,n):
    g_mn=np.sqrt((2*n+1)/(4*np.pi*n*(n+1)))*np.sqrt(factorial(n-m)/factorial(m+n))

    return g_mn

def bessel(x,n,kind=1,derivative=0):
    # upto 4th der 
    # for zero value, use non-zero instead
    x=x.detach().numpy()
    '''
    
        der_fuc=torch.sin(x)/x
        j_n_test=der_fuc
        for n_test in range(0,n):
            j_n_test=nth_derivative(j_n_test,x,1)
            j_n_test=j_n_test*(1/x)
        j_n=j_n_test*((-1)**n)*(x**n)

        
        '''
    if kind==1:
        J_n=spherical_jn(n, x,derivative=derivative)


        #j_n=nth_derivative(j_n,x,derivative)

    elif kind==2:
        J_n=spherical_yn(n, x,derivative=derivative)

    else: 
        pass

    return J_n

def hankel(x,n,kind=1,derivative=0):
    '''
    h_n_real=bessel(x,n,kind=1,tensor_form=True)
    h_n_imag=bessel(x,n,kind=2,tensor_form=True)
    h_n_real_der=nth_derivative(h_n_real, x, derivative)
    h_n_imag_der=nth_derivative(h_n_imag, x, derivative)
    h_n_real_der=h_n_real_der.detach().numpy()
    h_n_imag_der=h_n_imag_der.detach().numpy()
'''
    x=x.detach().numpy()
    if kind==1:
        h_n=bessel(x,n,kind=1,derivative=derivative)+bessel(x,n,kind=2,derivative=derivative)*1j
    elif kind==2:
        h_n=bessel(x,n,kind=1,derivative=derivative)-bessel(x,n,kind=2,derivative=derivative)*1j

    else: 
        pass    
    return h_n


def MN_mn(r,nu,phi,e_k,m,n,RG=True):
    k=np.sqrt(np.sum(e_k**2))
    kr=k*r
    g_mn=gamma_mn(m,n)
    CBP=CBP_mn(nu,phi,m,n,CBP_theta=True)
    C=CBP[0]
    B=CBP[1]
    P=CBP[2]
 



    if RG==True:
        radi_fuc=bessel(kr,n,kind=1,derivative=0)
 
        radi_fuc_der=bessel(kr,n,kind=1,derivative=1)
    else:
        radi_fuc=hankel(kr,n,kind=1,derivative=0)
        radi_fuc_der=hankel(kr,n,kind=1,derivative=1)
   
    radi_fuc=np.array([radi_fuc,radi_fuc,radi_fuc])
    
    radi_fuc_der=np.array([radi_fuc_der,radi_fuc_der,radi_fuc_der])

    
    M_mn=g_mn*C*radi_fuc


    # N =N_1+N_2

    kr=kr.detach().numpy()
    N_mn=g_mn*((n*(n+1)/kr)*radi_fuc*P+((1/kr)*radi_fuc+radi_fuc_der)*B)




    return M_mn,N_mn


def ab_matrix(n,m,e_polar,e_k):

    '''
    C.57,58

    a: MN=0
    b: MN=1
    '''
    angle_ek=angle_from_ek(e_k)
 
    theta_ek=torch.tensor(angle_ek[0],requires_grad=True)
    phi_ek=torch.tensor(angle_ek[1],requires_grad=True)
    CBP=CBP_mn(theta_ek,phi_ek,m,n,CBP_theta=False)
 
    
    C_star=np.conjugate(CBP[0])
    B_star=np.conjugate(CBP[1])


    d_n=np.sqrt((2*n+1)/(4*np.pi*n*(n+1)))

    phi_ek=phi_ek.detach().numpy()
    mn_factor=((1j)**(n-1))*4*np.pi*((-1)**m)*d_n*np.exp(-1j*m*phi_ek)

    a_mn=mn_factor*np.dot(e_polar,C_star)*(1j)
    b_mn=mn_factor*np.dot(e_polar,B_star)


    return a_mn,b_mn


def rec_plane_wave(r,nu,phi,e_k,e_polar,n_max):
    #e_k=e_k / (e_k**2).sum()**0.5

    size=r.shape
    E_total_x=np.zeros(size)
    E_total_y=np.zeros(size)
    E_total_z=np.zeros(size)

    
    for n in range(1,n_max):
 
        for m in range(-n,n+1):
            MN=MN_mn(r,nu,phi,e_k,m,n,RG=True)
            AB=ab_matrix(n,m,e_polar,e_k)
           


            M_x=MN[0][0,:,:]*AB[0]
            M_y=MN[0][1,:,:]*AB[0]
            M_z=MN[0][2,:,:]*AB[0]

            N_x=MN[1][0,:,:]*AB[1]
            N_y=MN[1][1,:,:]*AB[1]
            N_z=MN[1][2,:,:]*AB[1]

            E_total_x=E_total_x+M_x+N_x
            E_total_y=E_total_y+M_y+N_y
            E_total_z=E_total_z+M_z+N_z
 
    return E_total_x, E_total_y,E_total_z













#---------Plane Wave--------------------------

def plane_wave_dir(polar_psi,polar_chi,E0,k,direction=False):
    '''
    planewave with direction in k.
    x,y,z can be mesh and k is a vector. E0 id a scalar 
    '''
    # build local frame by using The Gram-Schmidt Process
    k_hat=k / (k**2).sum()**0.5
    # choose random dir 
    e_k=k_hat
    check_ort=e_k[np.abs(e_k)>0]

    
    if len(check_ort)>1:
        Index_ture=np.where(np.any(e_k!=0))
        e_a=np.array([1,1,1])
        bot_ea_2=e_k[Index_ture[0]]

        e_a_val=-(np.sum(e_k)-bot_ea_2)/bot_ea_2
        e_a[Index_ture[0]]=e_a_val

        e_a=e_a/ (e_a**2).sum()**0.5

        

    elif len(check_ort)==1:
        if np.abs(e_k[0])==1:
            e_a=np.array([0,0,1])
        elif np.abs(e_k[1])==1:
            e_a=np.array([0,0,-1])
        elif np.abs(e_k[2])==1:
            e_a=np.array([-1,0,0])
    # get third basis 
    e_b=np.cross(e_a,k_hat)
    #print(e_b)

 
    e_b=e_b/ (e_b**2).sum()**0.5

    # get the angle 
    E_a=E0*(np.cos(polar_chi)*np.sin(polar_psi)-1j*np.sin(polar_chi)*np.cos(polar_psi))
    E_b=E0*(np.cos(polar_chi)*np.cos(polar_psi)+1j*np.sin(polar_chi)*np.sin(polar_psi))


    e_polar=E_a*e_a+E_b*e_b

    if direction==True:
        return e_a,e_b
    else:

        return e_polar

def angle_from_ek(ek):
    ek=ek / (ek**2).sum()**0.5
 
    nu=np.arccos(ek[2])
    phi=np.arccos(ek[0]/(np.sin(nu)))
    '''
    if np.isnan(phi)==True:
        ek=[0,1,0]

    elif phi==0 and nu==0:
        ek=[0,0,1]
 
    else:      
        pass
 
    nu=np.arccos(ek[2])
    phi=np.arccos(ek[0]/(np.sin(nu)))
    '''
 
    return nu, phi

def plane_wave_grid(x,y,z,e_polar,K,omega=0,t0=0):
    '''
    generating grid for the inputing array
    '''
#    total_len_x=x.shape[0]*x.shape[1]*x.shape[2]
 #   total_len_y=y.shape[0]*y.shape[1]*y.shape[2]
  #  total_len_z=z.shape[0]*z.shape[1]*z.shape[2]
    k=np.sqrt(K[0]**2+K[1]**2+K[2]**2)
    wavelen=2*pi/k
 
    check_x_shape=x.shape
    check_y_shape=y.shape
    check_z_shape=z.shape

    x_r=x.reshape(-1,order='F')
    y_r=y.reshape(-1,order='F')
    z_r=z.reshape(-1,order='F')
    Total_field=np.array([x_r,y_r,z_r])
  
    #K_array=np.tile(K, (total_len_x,1))
    Dot_array=np.dot(Total_field.T,K)
 
    factor_exp=np.exp(1j*(Dot_array-omega*t0)).T
 
    
 
    E_field_x=e_polar[0]*factor_exp.reshape(check_x_shape,order='F')
    E_field_y=e_polar[1]*factor_exp.reshape(check_y_shape,order='F')
    E_field_z=e_polar[2]*factor_exp.reshape(check_z_shape,order='F')

    
    return E_field_x,E_field_y,E_field_z
    
     
def plane_sec_mesh(direction_vector1,direction_vector2,span1,span2,density):
    range_1=np.linspace(-span1,span1,density)
    range_2=np.linspace(-span2,span2,density)

    vector_1_x=range_1*direction_vector1[0]
    vector_1_y=range_1*direction_vector1[1]
    vector_1_z=range_1*direction_vector1[2]

    vector_2_x=range_2*direction_vector2[0]
    vector_2_y=range_2*direction_vector2[1]
    vector_2_z=range_2*direction_vector2[2]


    x_grid=np.meshgrid(vector_1_x,vector_2_x)
    y_grid=np.meshgrid(vector_1_y,vector_2_y)
    z_grid=np.meshgrid(vector_1_z,vector_2_z)

    return x_grid 






















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

    e_polar=np.array([E_a*e_a[0]+E_b*e_b[0],E_a*e_a[1]+E_b*e_b[1],E_a*e_a[2]+E_b*e_b[2]])

    return E_x,E_y,E_z,e_polar

# %%
