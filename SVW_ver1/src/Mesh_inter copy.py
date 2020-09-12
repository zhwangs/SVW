#%%
from SphTools import *
from GeometryTools import *
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
        Theta=np.linspace(0.00001,np.pi,N)
        Phi=np.linspace(0.00001,2*np.pi,N)
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

    def mesh_save_rec(self,save_path='/home/zhw/SphEM-1/test/ShoMesh_data/shpmesh.mat'):
        '''Save the file into matlab var '''
        Data_mesh = {
        "radius_mesh_rec": self.Radius_grid_NxN,
        "theta_mesh_rec": self.Theta_rec,
        "phi_mesh_rec": self.Phi_rec,
        "dir_mesh_rec": self.e_r_rec
        }
        sio.savemat(save_path, Data_mesh)



#%% Testing 
gmsh_file = '/home/zhw/SphEM-1/cache/sample_mesh/mesh_5.msh'
class BSVHE(Mesh_inter):
    '''
    Boundary Sphercial Vector Harmonic expansion (BSVHE). 
    The objective is to expand the electric field on a given geometric serface,
    in the basis of Vector Spherical Waves. 
    #############
    Variables:
    #interpo_density: the density of triangle interpolation: defalt=20,
    #mesh_grid_size: The target euqal space grid gmsh_file, Geo_size
    '''
    def __init__(self,gmsh_file,
                N_mesh_save_path='/home/zhw/SphEM-1/test/ShoMesh_data/shpmesh.mat',
                interpo_density=20,mesh_grid_size=20,
                        polar_psi=0,polar_chi=0,a=0,b=0,E0=1):
        super().__init__(gmsh_file)
        # interpolation (triangle mesh)
        Data_g=Mesh_inter.mesh_interpo(self,interpo_density)

        K=Mesh_inter.L_grid(self,mesh_grid_size)
        K_theta=K[3]
        K_phi=K[4]
        K_radius=K[2]
        R=K[1]
        K_dir=K[0]
        Rec=Mesh_inter.generate_rec_mesh(self,K_theta,K_phi,K_radius)
        Mesh_inter.mesh_save_rec(self,N_mesh_save_path)
        x_raw,y_raw,z_raw=Rec[0],Rec[1],Rec[2]
        self.x_raw=x_raw
        self.y_raw=y_raw
        self.z_raw=z_raw

        self.polar_psi=polar_psi
        self.polar_chi=polar_chi
        self.a=a
        self.b=b
        self.E0=E0

        


    def geo_reshape(self,Geo_size=3e-6):
        self.x=Geo_size*self.x_raw
        self.y=Geo_size*self.y_raw
        self.z=Geo_size*self.z_raw
        self.Geo_size=Geo_size


        

    
    def plane_wave_at_surf(self,polar_psi=0,polar_chi=0,a=0,b=0,E0=1,Omega=3e14):
        '''
        An example of a wave form at the geometry surface. 

        '''
        Mu_vacuum = 1.25663706212e-6    # vacuum magnetic permeability
        Eps_vacuum = 8.8541878128e-12   # vacuum electric permittivity
        ZVAC = 376.73031346177          # impedance of free space
        c = 299792458                   # the speed of light, m/s
        pi = np.pi
        Eps=Eps_vacuum
        Mu = Mu_vacuum

        k=Omega/c#2*pi/span 

        
        x=self.x
        y=self.y
        z=self.z
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

        self.E_x=E_x
        self.E_y=E_y
        self.E_z=E_z
        return E_x,E_y,E_z

    def surf_exp(self, E_x=0,E_y=0,E_z=0,L_max=2,Omega=3e14):

        Mu_vacuum = 1.25663706212e-6    # vacuum magnetic permeability
        Eps_vacuum = 8.8541878128e-12   # vacuum electric permittivity
        ZVAC = 376.73031346177          # impedance of free space
        c = 299792458                   # the speed of light, m/s
        pi = np.pi
        Eps=Eps_vacuum
        Mu = Mu_vacuum

        k=Omega/c#2*pi/span 

        totaL_len=len(self.x)**2
        x=self.x
        y=self.y
        z=self.z
        self.L_max=L_max
        r,theta,phi = convert_coordinates_Cart2Sph(x, y, z)
        theta=np.round(theta,5)
 

        N_max=L_max
        self.L_max=L_max
        Alpha_list=np.array([])
        # n loop 

        alpha_max=get_N_SphericalWaves_given_L_max(N_max)

        Row_vec=np.zeros(  (3*totaL_len, alpha_max+2)  )*1j  # I have just made 'alpha_max' to be at the column


        from progressbar import ProgressBar
        pbar = ProgressBar()
        for n in range(1,N_max+1):
            if n==1:
            # pass
                Row_vec[:,0]=np.ones(totaL_len*3)#+np.ones(totaL_len*3)*1j
                Row_vec[:,1]=np.ones(totaL_len*3)*1j
            #if n==2:
            #    Row_vec[:,1]=np.ones(totaL_len*3)*1j

            for m in range(-n,n+1):
                #print(m)
                #print(m)
                MN=get_M_N_waves(r,theta,phi, m,n, Omega, Eps, Mu, kind=1, return_Cartesian=True)

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

        if (E_x+E_y+E_z)==0: 
            #print('Plane-Wave Demo')
            E_x_s=self.E_x
            E_y_s=self.E_y
            E_z_s=self.E_z
        else:
            E_x_s=E_x
            E_y_s=E_y
            E_z_s=E_z           

        
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

        self.leng_cal=leng_cal
        self.E_cal_real_real=E_cal_real_real
        self.E_p_real=E_p_real
        self.E_cal_imag_real=E_cal_imag_real
        self.E_p_imag=E_p_imag
        self.E_p_x_r=E_p_x_r
        self.E_p_y_r=E_p_y_r
        self.E_p_z_r=E_p_z_r
        self.E_x_r=E_x_r
        self.E_y_r=E_y_r
        self.E_z_r=E_z_r
        self.E_mag_r=E_mag_r


        # Calculate Error 
        E_real_ave_error_precent=np.sqrt(np.mean((E_cal_real_real-E_p_real)**2))/np.sqrt(np.mean((E_cal_real_real+E_p_real)**2))
        E_imag_ave_error_precent=np.sqrt(np.mean((E_cal_imag_real-E_p_imag)**2))/np.sqrt(np.mean((E_cal_imag_real+E_p_imag)**2))

        E_real_span_ratio=np.sign(np.max(E_cal_real_real)-np.max(E_p_real))*(((np.max(E_cal_real_real)-np.min(E_cal_real_real))-(np.max(E_p_real)+np.min(E_p_real)))/((np.max(E_cal_real_real)-np.min(E_cal_real_real))+(np.max(E_p_real)-np.min(E_p_real))))
        E_imag_span_ratio=np.sign(np.max(E_cal_imag_real)-np.max(E_p_imag))*(((np.max(E_cal_imag_real)-np.min(E_cal_imag_real))-(np.max(E_p_imag)+np.min(E_p_imag))))/((np.max(E_cal_imag_real)-np.min(E_cal_imag_real))+(np.max(E_p_imag)-np.min(E_p_imag)))
        
        return E_real_ave_error_precent,E_imag_ave_error_precent,E_real_span_ratio,E_imag_span_ratio

    
    def see_exp(self, index_plot=False,p_plot=False,d_plot=False):
        import numpy as np
        from matplotlib import pyplot as plt
        import matplotlib.cm as cm
        from matplotlib.colors import Normalize
        from mpl_toolkits.mplot3d import Axes3D

        if index_plot:
            f,axs=plt.subplots(2,3)
            place_hod=np.linspace(0,1,int(self.leng_cal/2))
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
            f.suptitle('Plane-Wave reco by index(geo-size='+str(self.Geo_size)+'), L_max='+ str(self.L_max), fontsize=30)

            ax_1.plot(place_hod,self.E_cal_real_real)
            ax_1.set_title('Calculated Real')
            ax_1.set_xlabel('normalized index ')
            ax_1.set_ylabel('Field Strength (E.real)')
            ax_2.plot(place_hod,self.E_p_real)
            ax_2.set_title('Plane-Wave Real')
            ax_2.set_xlabel('normalized index ')
            ax_2.set_ylabel('Field Strength (E.real)')
            ax_6.plot(place_hod,np.abs(self.E_cal_real_real-self.E_p_real)**2)
            ax_6.set_title('Difference,Real')
            ax_6.set_xlabel('normalized index ')
            ax_6.set_ylabel('Field Strength Difference (E_rec.real-E_plane.real)')

            ax_3.plot(place_hod,(self.E_cal_imag_real))
            ax_3.set_title('Calculated imag')
            ax_3.set_xlabel('normalized index ')
            ax_3.set_ylabel('Field Strength (E.imag)')
            ax_4.plot(place_hod,self.E_p_imag)
            ax_4.set_title('Plane-Wave imag')
            ax_4.set_xlabel('normalized index ')
            ax_4.set_ylabel('Field Strength (E.imag)')
            ax_5.plot(place_hod,(np.abs(self.E_cal_imag_real-self.E_p_imag)**2))
            ax_5.set_title('Difference, imag')
            ax_5.set_xlabel('normalized index ')
            ax_5.set_ylabel('Field Strength Difference (E_rec.imag-E_plane.imag)')

        if p_plot:
            fig = plt.figure()
            fig.suptitle('Plane-Wave reco (geo-size='+str(self.Geo_size)+'), L_max='+ str(self.L_max), fontsize=30)
            ax_1 = fig.add_subplot(221)
            ax_2 = fig.add_subplot(222)
            ax_3 = fig.add_subplot(223)
            ax_4 = fig.add_subplot(224, projection='3d')


            fig.set_figheight(15)
            fig.set_figwidth(15)

            ax_1.quiver(self.x,self.y,self.E_x_r ,self.E_y_r)
            ax_1.set_title('x-y plane')
            ax_1.set_xlabel('x (meter)')
            ax_1.set_ylabel('y (meter)')

            ax_2.quiver(self.x,self.z,self.E_x_r ,self.E_z_r)
            ax_2.set_title('x-z plane')
            ax_2.set_xlabel('x (meter)')
            ax_2.set_ylabel('z (meter)')

            ax_3.quiver(self.y,self.z,self.E_y_r ,self.E_z_r)
            ax_3.set_title('y-z plane')
            ax_3.set_xlabel('y (meter)')
            ax_3.set_ylabel('z (meter)')




            norm = Normalize()
            colors = norm(self.E_mag_r**2) # auto-adjust true radius into [0,1] for color mapping

            cmap = cm.get_cmap("coolwarm")
            ax_4.plot_surface(self.x, self.y,self.z, linewidth=0, facecolors=cmap(colors), shade=True, alpha=1)

            # the surface is not mappable, we need to handle the colorbar manually
            mappable = cm.ScalarMappable(cmap=cmap)
            mappable.set_array(self.E_mag_r**2)
            fig.colorbar(mappable, shrink=0.2, aspect=5,label='Field strength (E^2)')
            # auto-adjust true radius into [0,1] for color mapping
            ax_4.quiver(self.x, self.y, self.z, self.E_x_r, self.E_y_r, self.E_z_r, length=0.2*self.x.max())
            ax_4.view_init(elev=15, azim=45)
            ax_4.set_title('3D view, vector field')
            ax_4.set_xlabel('x (meter)')
            ax_4.set_ylabel('y (meter)')
            ax_4.set_zlabel('z (meter)')

        if d_plot:
            fig = plt.figure()
            fig.set_figheight(15)
            fig.set_figwidth(15)

            ax = fig.add_subplot(111, projection='3d')
            ax.view_init(elev=15, azim=45)#np.pi/4)

            norm = Normalize()
            colors = norm(self.E_mag_r) # auto-adjust true radius into [0,1] for color mapping

            cmap = cm.get_cmap("coolwarm")
            ax.plot_surface(self.x, self.y, self.z, linewidth=0, facecolors=cmap(colors), shade=True, alpha=1)

            # the surface is not mappable, we need to handle the colorbar manually
            mappable = cm.ScalarMappable(cmap=cmap)
            mappable.set_array(self.E_mag_r**2)
            fig.colorbar(mappable, shrink=0.2, aspect=5,label='Field strength (E^2)')
            # auto-adjust true radius into [0,1] for color mapping
        
            ax.quiver(self.x, self.y, self.z, self.E_x_r, self.E_y_r, self.E_z_r, length=0.2*self.Geo_size)

            ax.set_title('3D view, vector field')
            ax.set_xlabel('x (meter)')
            ax.set_ylabel('y (meter)')
            ax.set_zlabel('z (meter)')


    def geo_freq_error_map(self, wave_len_start,wave_len_end,Geo_size_start,Geo_size_end,mesh_density):
        # error map that gives error respect to different frequency and geometric size. 
        from progressbar import ProgressBar
        pbar = ProgressBar()
        wave_len=np.linspace(wave_len_start,wave_len_end,mesh_density)
        # frequency
        c = 299792458 
        omega=2*np.pi*c/wave_len
        geo_size=np.linspace(Geo_size_start,Geo_size_end,mesh_density)
        wave_len_grid,geo_grid=np.meshgrid(wave_len,geo_size)
        Error_grad_real=np.zeros((mesh_density,mesh_density))
        Error_grad_imag=np.zeros((mesh_density,mesh_density))
        Error_grad_real_bias=np.zeros((mesh_density,mesh_density))
        Error_grad_imag_bias=np.zeros((mesh_density,mesh_density))

        F=BSVHE(gmsh_file,interpo_density=20,mesh_grid_size=20)
        for i in pbar(range(0,mesh_density)):
            current_omega=omega[i]
            current_geo_size=geo_size[i]

            for j in range(0,mesh_density):
                S=F.plane_wave_at_surf(polar_psi=self.polar_psi,polar_chi=self.polar_chi,a=self.a,b=self.b,E0=self.E0,Omega=omega[j])
           
                F.geo_reshape(Geo_size=current_geo_size)
                K=F.surf_exp(Omega=current_omega)
                Error_grad_real[i,j]=K[0]
                Error_grad_imag[i,j]=K[1]
                Error_grad_real_bias[i,j]=K[2]
                Error_grad_imag_bias[i,j]=K[3]

        self.wave_len_grid= (wave_len_grid)
        self.geo_grid= (geo_grid)
        self.Error_grad_real=Error_grad_real
        self.Error_grad_imag=Error_grad_imag
        self.Error_grad_real_bias=Error_grad_real_bias
        self.Error_grad_imag_bias=Error_grad_imag_bias
        self.wave_len_start=wave_len_start
        self.wave_len_end=wave_len_end
        self.Geo_size_start=Geo_size_start
        self.Geo_size_end=Geo_size_end
        self.mesh_density=mesh_density
        return wave_len_grid,geo_grid,Error_grad_real,Error_grad_imag,Error_grad_real_bias,Error_grad_imag_bias


    def see_error_map(self,geo_size_show=False,geo_step_show=False):
        if geo_size_show:
            wave_len_grid=self.wave_len_grid
            geo_size=self.geo_grid
            Error_grad_real=self.Error_grad_real
            Error_grad_imag=self.Error_grad_imag
            Error_grad_real_bias=self.Error_grad_real_bias
            Error_grad_imag_bias=self.Error_grad_imag_bias
            
            fig = plt.figure()
            fig.suptitle('Error map', fontsize=30)
            ax_1 = fig.add_subplot(221)
            ax_2 = fig.add_subplot(222)
            ax_3 = fig.add_subplot(223)
            ax_4 = fig.add_subplot(224)

            fig.set_figheight(15)
            fig.set_figwidth(15)

            a1=ax_1.pcolormesh(wave_len_grid, geo_size,Error_grad_real)
            fig.colorbar(a1, ax=ax_1)
            ax_1.set_title('Real Error heat map')
            ax_1.set_xlabel('wave length (m)')
            ax_1.set_ylabel(' geometric size (m)')
            a2=ax_2.pcolormesh(wave_len_grid, geo_size,Error_grad_imag)
            fig.colorbar(a2, ax=ax_2)
            ax_2.set_title('Imag Error heat map')
            ax_2.set_xlabel('wave length (m)')
            ax_2.set_ylabel(' geometric size (m)')

            a3=ax_3.pcolormesh(wave_len_grid, geo_size,Error_grad_real_bias)
            fig.colorbar(a3, ax=ax_3)
            ax_3.set_title('Real Error bias heat map')
            ax_3.set_xlabel('wave length (m)')
            ax_3.set_ylabel(' geometric size (m)')
            a4=ax_4.pcolormesh(wave_len_grid, geo_size,Error_grad_imag_bias)
            fig.colorbar(a4, ax=ax_4)
            ax_4.set_title('Imag Error bias heat map')
            ax_4.set_xlabel('wave length (m)')
            ax_4.set_ylabel(' geometric size (m)')
            #fig.text( 0,0," Wave length range: ["+str(self.wave_len_start)+','+str(self.wave_len_end)+']'
          #          +" Geo-size range: ["+str(self.Geo_size_start)+','+str(self.Geo_size_end)+']'+', grid density='+str(self.mesh_density), size=30)
        if geo_step_show:
            wave_len_grid=self.precent_freq_grid
            geo_size=self.geo_grid_1
            Error_grad_real=self.Error_grad_real_1
            Error_grad_imag=self.Error_grad_imag_1
            Error_grad_real_bias=self.Error_grad_real_bias_1
            Error_grad_imag_bias=self.Error_grad_imag_bias_1
            

            fig = plt.figure()
            fig.suptitle('Error map', fontsize=30)
            ax_1 = fig.add_subplot(221)
            ax_2 = fig.add_subplot(222)
            ax_3 = fig.add_subplot(223)
            ax_4 = fig.add_subplot(224)

            fig.set_figheight(15)
            fig.set_figwidth(15)

            a1=ax_1.pcolormesh(wave_len_grid, geo_size,Error_grad_real)
            fig.colorbar(a1, ax=ax_1)
            ax_1.set_title('Real Error heat map')
            ax_1.set_xlabel('wave length (m)')
            ax_1.set_ylabel(' geometric size (m)')
            a2=ax_2.pcolormesh(wave_len_grid, geo_size,Error_grad_imag)
            fig.colorbar(a2, ax=ax_2)
            ax_2.set_title('Imag Error heat map')
            ax_2.set_xlabel('wave length (m)')
            ax_2.set_ylabel(' geometric size (m)')

            a3=ax_3.pcolormesh(wave_len_grid, geo_size,Error_grad_real_bias)
            fig.colorbar(a3, ax=ax_3)
            ax_3.set_title('Real Error bias heat map')
            ax_3.set_xlabel('wave length (m)')
            ax_3.set_ylabel(' geometric size (m)')
            a4=ax_4.pcolormesh(wave_len_grid, geo_size,Error_grad_imag_bias)
            fig.colorbar(a4, ax=ax_4)
            ax_4.set_title('Imag Error bias heat map')
            ax_4.set_xlabel('wave length (m)')
            ax_4.set_ylabel(' geometric size (m)')
            #fig.text( 0,0," Wave length range: ["+str(self.wave_len_start)+','+str(self.wave_len_end)+']'
             #       +" Geo-size range: ["+str(self.Geo_size_start)+','+str(self.Geo_size_end)+']'+', grid density='+str(self.mesh_density), size=30)            
    
    def geo_step_error_map(self,Geo_size_start,Geo_size_end,mesh_density):
        # error map that gives error respect to different frequency and geometric size. 
        from progressbar import ProgressBar
        pbar = ProgressBar()
        

        geo_size=np.linspace(Geo_size_start,Geo_size_end,mesh_density)
        precent_freq_grid,geo_grid_1=np.meshgrid(geo_size,geo_size)
        c = 299792458 
        geo_eq_freq=2*np.pi*c/geo_size
        #geo_grid_1=np.tile(geo_eq_freq,(mesh_density,1)).T
        
        #precent_freq_grid,geo_grid=np.meshgrid(precent_freq,geo_size)
        Error_grad_real_1=np.zeros((mesh_density,mesh_density))
        Error_grad_imag_1=np.zeros((mesh_density,mesh_density))
        Error_grad_real_bias_1=np.zeros((mesh_density,mesh_density))
        Error_grad_imag_bias_1=np.zeros((mesh_density,mesh_density))

        #Precent_freq_grid=np.zeros((mesh_density,mesh_density))

        F=BSVHE(gmsh_file,interpo_density=20,mesh_grid_size=20)
        for i in pbar(range(0,mesh_density)):
            current_geo_size=geo_size[i]
            #precent_freq=current_geo_size*np.linspace(0.1,100,mesh_density)
            # frequency
            #Precent_freq_grid[i,:]=precent_freq
            omega=geo_eq_freq
            
            F.geo_reshape(Geo_size=current_geo_size)

            S=F.plane_wave_at_surf(polar_psi=self.polar_psi,polar_chi=self.polar_chi,a=self.a,b=self.b,E0=self.E0, Omega=geo_eq_freq[i])
           
            for j in range(0,mesh_density):
                current_omega=omega[j]

                K=F.surf_exp(Omega=current_omega)
                Error_grad_real_1[i,j]=K[0]
                Error_grad_imag_1[i,j]=K[1]
                Error_grad_real_bias_1[i,j]=K[2]
                Error_grad_imag_bias_1[i,j]=K[3]

        self.precent_freq_grid= (precent_freq_grid)
        self.geo_grid_1= (geo_grid_1)
        self.Error_grad_real_1=Error_grad_real_1
        self.Error_grad_imag_1=Error_grad_imag_1
        self.Error_grad_real_bias_1=Error_grad_real_bias_1
        self.Error_grad_imag_bias_1=Error_grad_imag_bias_1


        self.Geo_size_start=Geo_size_start
        self.Geo_size_end=Geo_size_end
        self.mesh_density=mesh_density
        return precent_freq_grid,geo_grid_1,Error_grad_real_1,Error_grad_imag_1,Error_grad_real_bias_1,Error_grad_imag_bias_1

#%%

F=BSVHE(gmsh_file,interpo_density=20,mesh_grid_size=20,polar_psi=1,polar_chi=1,a=0,b=1,E0=1)
#F.geo_reshape(Geo_size=3e-6)
#S=F.plane_wave_at_surf(0,0,0,0,1)
wave_len_start=1e-6
wave_len_end=1e-8
Geo_size_start=1e-6
Geo_size_end=1e-8

mesh_density=30

Error=F.geo_freq_error_map(wave_len_start,wave_len_end,Geo_size_start,Geo_size_end,mesh_density)
F.see_error_map(geo_size_show=True)

#%% 
Geo_map=F.geo_step_error_map(Geo_size_start,Geo_size_end,mesh_density)
#%%
# Geo_map[1]
F.see_error_map(geo_step_show=True)




#%%
polar_psi=1
polar_chi=1
a=0
b=1
E0=1
F=BSVHE(gmsh_file,interpo_density=20,mesh_grid_size=20,polar_psi=polar_psi,polar_chi=polar_chi,a=a,b=b,E0=E0)

F.geo_reshape(Geo_size=3e-6)
S=F.plane_wave_at_surf(polar_chi=polar_chi,polar_psi=polar_psi,a=a,b=b,E0=E0)


K=F.surf_exp()
F.see_exp(index_plot=True,p_plot=True,d_plot=True)


# %%


# %%


# %%
fig = plt.figure()
fig.suptitle('Error map', fontsize=30)


# %%
