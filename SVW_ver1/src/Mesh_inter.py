#%%
#from SphTools import *
#from GeometryTools import *
import numpy as np
import scipy as sp
 
import sys
import pandas as pd

from mpl_toolkits import mplot3d
import scipy.io as sio
#import pyshtools
class Mesh_inter:
    def __init__(self,gmsh_file):
        import meshio 
        # mesh file 
        mesh_object=meshio.read(gmsh_file)
        mesh_points=mesh_object.points
        mesh_triangle=mesh_object.cells_dict['triangle']
        # recenter the mesh to its geo-center 
        mesh_points_center=np.mean(mesh_points,axis=0)
        mesh_points=mesh_points-mesh_points_center
        # access triangle length
        mesh_len=len(mesh_triangle[:,0])
            # Reshape the mesh in to a 1D array
        r=mesh_triangle.reshape((1,3*mesh_len))
        mesh_data_list=np.array(mesh_points[r[0]])
        mesh_data=mesh_data_list.reshape((mesh_len,3,3))
        # center location for triangles 
        mesh_tri_center=np.mean(mesh_data,axis=1)
        # Reconstract matrix to 1D array
        mesh_data_tri=mesh_data_list.reshape((int(mesh_len),9))

        self.mesh_data_tri=mesh_data_tri
        self.mesh_tri_center=mesh_tri_center
        self.mesh_data=mesh_data
        self.mesh_len=mesh_len
        self.mesh_points=mesh_points
       
        

    def mesh_interpo(self,mesh_density):
        ''' A function interpolate triangular mesh with density N=mesh_density'''
        # parameter equation
        self.mesh_density=mesh_density
        u=np.linspace(0,1,mesh_density,endpoint=False)
        # generate mesh for v (each row represents v: 0 to 1-u)
        end_v=np.ones(len(u))
        start_v=np.zeros(len(u))
        v=np.linspace(start_v,end_v-u,mesh_density,endpoint=False).T
        # Stack and repeat same v for mesh_len (number of triangles) times
        v=np.tile(v,(self.mesh_len,1))
        # generate mesh for u (each row represent one u)
        u=np.repeat(np.array([u]),mesh_density, axis=0).T
        u=np.tile(u,(self.mesh_len,1))

        
        def tri_para_fuc(u,v):
            ''' A function interpolate triangular mesh with density N=mesh_density
            # Constract matrix for triangle parametric equations'''
            # Coff matrix (9,1) A+Bu+Cv
            A_x=np.array([1,0,0,0,0,0,0,0,0])
            A_y=np.array([0,1,0,0,0,0,0,0,0])
            A_z=np.array([0,0,1,0,0,0,0,0,0])

            B_x=np.array([-1,0,0,1,0,0,0,0,0])
            B_y=np.array([0,-1,0,0,1,0,0,0,0])
            B_z=np.array([0,0,-1,0,0,1,0,0,0])
            
            C_x=np.array([-1,0,0,0,0,0,1,0,0])
            C_y=np.array([0,-1,0,0,0,0,0,1,0])
            C_z=np.array([0,0,-1,0,0,0,0,0,1])
            #check data uniqueness
            #Difference_check=self.mesh_data_tri[:,1]-self.mesh_data_tri[:,0]-self.mesh_data_tri[:,2]
            #Unique_diff=np.unique(Difference_check)

            # Shape transform 
            A_dot_x=np.repeat(np.array([np.dot(self.mesh_data_tri,A_x.T)]),mesh_density, axis=0).T
            A_dot_x=np.repeat(A_dot_x,mesh_density, axis=0)
            B_dot_x=np.repeat(np.array([np.dot(self.mesh_data_tri,B_x.T)]),mesh_density, axis=0).T
            B_dot_x=np.repeat(B_dot_x,mesh_density, axis=0)
            C_dot_x=np.repeat(np.array([np.dot(self.mesh_data_tri,C_x.T)]),mesh_density, axis=0).T
            C_dot_x=np.repeat(C_dot_x,mesh_density, axis=0)

            A_dot_y=np.repeat(np.array([np.dot(self.mesh_data_tri,A_y.T)]),mesh_density, axis=0).T
            A_dot_y=np.repeat(A_dot_y,mesh_density, axis=0)
            B_dot_y=np.repeat(np.array([np.dot(self.mesh_data_tri,B_y.T)]),mesh_density, axis=0).T
            B_dot_y=np.repeat(B_dot_y,mesh_density, axis=0)
            C_dot_y=np.repeat(np.array([np.dot(self.mesh_data_tri,C_y.T)]),mesh_density, axis=0).T
            C_dot_y=np.repeat(C_dot_y,mesh_density, axis=0)

            A_dot_z=np.repeat(np.array([np.dot(self.mesh_data_tri,A_z.T)]),mesh_density, axis=0).T
            A_dot_z=np.repeat(A_dot_z,mesh_density, axis=0)
            B_dot_z=np.repeat(np.array([np.dot(self.mesh_data_tri,B_z.T)]),mesh_density, axis=0).T
            B_dot_z=np.repeat(B_dot_z,mesh_density, axis=0)
            C_dot_z=np.repeat(np.array([np.dot(self.mesh_data_tri,C_z.T)]),mesh_density, axis=0).T
            C_dot_z=np.repeat(C_dot_z,mesh_density, axis=0)

            #Repeat the calculation for different matrix order 
           
            # u vs v
            Para_x_uv=A_dot_x+B_dot_x*u+C_dot_x*v
            Para_y_uv=A_dot_y+B_dot_y*u+C_dot_y*v
            Para_z_uv=A_dot_z+B_dot_z*u+C_dot_z*v

            # v vs u 
            Para_x_vu=A_dot_x+B_dot_x*v+C_dot_x*u
            Para_y_vu=A_dot_y+B_dot_y*v+C_dot_y*u
            Para_z_vu=A_dot_z+B_dot_z*v+C_dot_z*u
            # append two calculations
            Para_x=np.append(Para_x_vu,Para_x_uv,axis=0)
            Para_y=np.append(Para_y_vu,Para_y_uv,axis=0)
            Para_z=np.append(Para_z_vu,Para_z_uv,axis=0)
 
            self.Para_x=Para_x
            self.Para_y=Para_y
            self.Para_z=Para_z
            self.Data_grid=np.array([Para_x,Para_y,Para_z])
            # return Mesh grid
            return Para_x,Para_y,Para_z
        # Call function
        Data_grid=tri_para_fuc(u,v)

        return Data_grid[0],Data_grid[1],Data_grid[2],self.mesh_points



    def angle_radius_mesh(self):

        ''' coordinate transformation: (x,y,z)->(r,theta,phi)
        Updated: nolonger needed for the process'''

        radius_mesh= np.sqrt(self.Para_x**2+ self.Para_y**2+self.Para_z**2)
        theta_mesh=np.arccos(self.Para_z/radius_mesh)
        phi_mesh=2*np.arctan(self.Para_y/self.Para_x)+np.pi

        self.radius_mesh=radius_mesh
        self.theta_mesh=theta_mesh
        self.phi_mesh=phi_mesh

        return radius_mesh,theta_mesh,phi_mesh

    def unit_direction_mesh(self):
        ''' coordinate unit vector reconstruction: (r,theta,phi)->(x,y,z)
        Updated: nolonger needed for the process'''
        e_r=np.array([np.sin(self.theta_mesh)*np.cos(self.phi_mesh),np.sin(self.theta_mesh)*np.sin(self.phi_mesh),np.cos(self.theta_mesh)])
 
        self.e_r=e_r
        return self.radius_mesh*e_r

    def mesh_save(self):
        ''' Save polar mesh file
        Updated: nolonger needed for the process '''     
        Data_mesh = {
        "radius_mesh": self.radius_mesh,
        "theta_mesh": self.theta_mesh,
        "phi_mesh": self.phi_mesh,
        "dir_mesh": self.e_r
        }
        sio.savemat('test.mat', Data_mesh)


    def G_quad_mesh_N_esti(self):
       '''Checking and validating the theta/phi increment Not used in mesh generation'''
       mesh_len=self.mesh_len
       mesh_density=self.mesh_density
       theta_mesh=self.theta_mesh
       phi_mesh=self.phi_mesh
       # reshape for theta 
       list_wise_theta=np.reshape(theta_mesh,2*mesh_len*mesh_density**2) # 3 (x,y,z), 2( u-v, v-u)
       # reshape for phi
       list_wise_phi=np.reshape(phi_mesh,2*mesh_len*mesh_density**2)
       # sort theta from 0-pi
       sorted_theta_list=np.sort(list_wise_theta)
       # record the index for theta sort
       sorted_theta_index_list=np.argsort(list_wise_theta)
       # apply to phi which find the coresponding phi
       sorted_phi_list_resp_theta=list_wise_phi[sorted_theta_index_list]

       # find unique array
       unique_theta_list=np.unique(sorted_theta_list, return_index=True,return_counts=True)
       unique_theta_arrray=unique_theta_list[0]
       unique_theta_index=unique_theta_list[1]
       unique_theta_counts=unique_theta_list[2]
       shifted_unique_theta_arrray=np.roll(unique_theta_arrray, -1)
       delta_theta=shifted_unique_theta_arrray-unique_theta_arrray
       s=np.sum(delta_theta[0:len(delta_theta)-1])
       # mean estimate datapoints 
       k_theta=np.pi/np.mean(delta_theta[0:len(delta_theta)-1])
       
       # Phi
       sorted_phi_list=np.sort(list_wise_phi)
       # record the index for phi sort
       sorted_phi_index_list=np.argsort(list_wise_phi)

       # apply to theta which find the coresponding theta
       sorted_theta_list_resp_theta=list_wise_theta[sorted_phi_index_list]
    
       # find unique array
       unique_phi_list=np.unique(sorted_phi_list, return_index=True,return_counts=True)
       unique_phi_arrray=unique_phi_list[0]
       unique_phi_index=unique_phi_list[1]
       unique_phi_counts=unique_phi_list[2]
       shifted_unique_phi_arrray=np.roll(unique_phi_arrray, -1)
       delta_phi=shifted_unique_phi_arrray-unique_phi_arrray
       s=np.sum(delta_phi[0:len(delta_phi)-1])
       # mean estimate datapoints 
       k_phi=2*np.pi/np.mean(delta_phi[0:len(delta_phi)-1])


       return delta_theta,k_theta,delta_phi,k_phi


    def GL_mesh(self, N):
        '''Looping attempts (slow): Updated: nolonger needed'''
        # checking the order of mag
        order_of_mag=int(np.floor(np.log10(N)))

        # equal space grid
        Theta=np.linspace(0,np.pi,N)
        Phi=np.linspace(0,2*np.pi,N)
        Theta=np.round(Theta,order_of_mag)
        Phi=np.round(Phi,order_of_mag)
        
        # range deduct
        theta_mesh=np.round(self.theta_mesh,order_of_mag)
        phi_mesh=np.round(self.phi_mesh,order_of_mag)
        radius_mesh=self.radius_mesh

        radius_array_targeted=[]
        for i in range(0,N):
            # grab theta from array
            theta_value_i=Theta[i]
            #Checking location in mash
            Theta_mesh_check=theta_mesh==theta_value_i
            # finding the radius and phi mesh that corresponding to the theta
            radius_mesh_theta_i=radius_mesh[Theta_mesh_check]
            phi_mesh_theta_i=phi_mesh[Theta_mesh_check]
            #print(radius_mesh_theta_i.shape)
            #print(phi_mesh_theta_i)
            for j in range(0,N):
                phi_value_j=Phi[j]
                Phi_submesh_check=phi_mesh_theta_i==phi_value_j
                # find the value for r for given phi_j and theta_i
                radius_submesh_phi_j=np.mean(radius_mesh_theta_i[Phi_submesh_check])
                #print(radius_submesh_phi_j)
                radius_array_targeted.append(radius_submesh_phi_j)

        Theta,Phi=np.meshgrid(Theta,Phi)
        radius_array_targeted=np.array(radius_array_targeted)
        radius_array_targeted=radius_array_targeted.reshape((N,N))

        return Theta,Phi,radius_array_targeted
                

    def generate_rec_mesh(self,Theta,Phi,R):
        '''reconstruction (x,y,z) from (Theta,Phi,R) '''
        e_r=np.array([np.sin(Theta)*np.cos(Phi),np.sin(Theta)*np.sin(Phi),np.cos(Theta)])
        #print(np.sin(self.theta_mesh))
        return R*e_r     


    def G_quad_mesh_N(self):
        ''' Second attempt without looping. Updated: function nolonger needed '''

       # N is the quad_mesh density
       # # define the mesh density for reshaping 
        mesh_len=self.mesh_len
        mesh_density=self.mesh_density
        theta_mesh=self.theta_mesh
        phi_mesh=self.phi_mesh
        radius_mesh=self.radius_mesh
        # reshape for theta 
        list_wise_theta=np.reshape(theta_mesh,2*mesh_len*mesh_density**2) # 3 (x,y,z), 2( u-v, v-u)
        # reshape for phi
        list_wise_phi=np.reshape(phi_mesh,2*mesh_len*mesh_density**2) 
        # reshape radius
        list_wise_radius=np.reshape(radius_mesh,2*mesh_len*mesh_density**2) 
        # sort theta from 0-pi
        sorted_theta_list=np.sort(list_wise_theta)
        # record the index for theta sort
        sorted_theta_index_list=np.argsort(list_wise_theta)

        # apply to phi and radius which find respects to theta 
        sorted_phi_list_resp_theta=list_wise_phi[sorted_theta_index_list]
        sorted_radius_list_resp_theta=list_wise_radius[sorted_theta_index_list]

        # reshape back
        theta_sorted_updated=np.reshape(sorted_theta_list,(2*mesh_len*mesh_density,mesh_density))
        phi_sorted_updated=np.reshape(sorted_phi_list_resp_theta,(2*mesh_len*mesh_density,mesh_density))
        radius_sorted_updated=np.reshape(sorted_radius_list_resp_theta,(2*mesh_len*mesh_density,mesh_density))
        return theta_sorted_updated,phi_sorted_updated,radius_sorted_updated
            


    def L_grid(self,N):
        ''' Equally spaced NxN (grid) construction'''
        Data_grid=self.Data_grid
        mesh_len=self.mesh_len
        mesh_density=self.mesh_density
        # create mesh for theta and phi 
        Theta=np.linspace(0.001,0.999*np.pi,N)
        Phi=np.linspace(0.001,1.999*np.pi,N)
        # generate mesh for Theta (each row represent one Theta)
        Theta=np.repeat(np.array([Theta]),N, axis=0).T

        # for phi
        Phi=np.repeat(np.array([Phi]),N, axis=0)
        # calculate the unit direction based on the theta and phi
        e_r=np.array([np.sin(Theta)*np.cos(Phi),np.sin(Theta)*np.sin(Phi),np.cos(Theta)])
        # reshape the unit direction into a (3,N^2) array for later dot product
        e_r=np.reshape(e_r,(3,N**2))
        # calcuate the radius mesh from the data grid found above
        Radius_grid=np.sqrt(Data_grid[0]**2+Data_grid[1]**2+Data_grid[2]**2)
        # normalize the data array
        Data_grid_normal=Data_grid/Radius_grid
        # reshape the data array into the form (3,2*mesh_len*mesh_density**2)
        Data_grid_normal=np.reshape(Data_grid_normal,(3,2*mesh_len*mesh_density**2))
        Data_grid=np.reshape(Data_grid,(3,2*mesh_len*mesh_density**2))
        print(e_r.shape)
        print(Data_grid_normal.shape)
        #Calculate the projection matrix for given data sets and targets. find min, so used 1-dot
        Projection_matrix=1-np.dot(e_r.T,Data_grid_normal)
        print(Projection_matrix.shape)
        mini_align_value=np.amin(Projection_matrix, axis=1)
        mini_align_arg_value=np.argmin(Projection_matrix, axis=1)
        # Find the most aligned direction 
        Data_grid_align=Data_grid.T[mini_align_arg_value]
        # Reshape the grid to NxN
        Radius_grid_NxN=np.sqrt(Data_grid_align[:,0]**2+Data_grid_align[:,1]**2+Data_grid_align[:,2]**2)
        Radius_grid_NxN=np.reshape(Radius_grid_NxN,(N,N))
        # Uniform N by N mesh
        print(Radius_grid_NxN.shape)

        self.Radius_grid_NxN=Radius_grid_NxN
        self.e_r_rec=e_r
        self.Theta_rec=Theta
        self.Phi_rec=Phi

        return e_r, Data_grid_align,Radius_grid_NxN,Theta,Phi

    def mesh_save_rec(self,path):
        '''Save the file into matlab var '''
        Data_mesh = {
        "radius_mesh_rec": self.Radius_grid_NxN,
        "theta_mesh_rec": self.Theta_rec,
        "phi_mesh_rec": self.Phi_rec,
        "dir_mesh_rec": self.e_r_rec
        }
        sio.savemat(path, Data_mesh)

