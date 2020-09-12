'''
We will pile all useful functions for handling geometry here
'''

# importing all functions from SphTools
from SphTools import *

# tool to read the mesh
import meshio # https://pypi.org/project/meshio/

# tool to convert map of r as a function of theta and pho into spherical harmonics
#import pyshtools



def get_goemetric_center(*args):
    '''
    Calculate the geometric center of a cloud of points
    
    Parameters:
    ===========
    you can either input one of the following:
        x,y,z       : ndarrays, represent the x,y,z points coordinates in Cartesian coordinates
    or
        points_array : ndarray, the last exis represent the x,y,z with indices 0,1,2, respectively.
    returns:
    ========
    geometric_center : ndarray, represent the geometric center. First, second and last element represent x,y and z coordinates.
    '''

    # to make our code flexible and accept different representations of point
    if len(args) == 1:
        x = args[0][...,0]
        y = args[0][...,1]
        z = args[0][...,2]

    elif len(args) == 3:
        x = args[0]
        y = args[1]
        z = args[2]
    
    else:
        print("there should be error here")
        # TODO: code an error here

    # calculation of geometric center
    geometric_center = np.zeros(3)
    geometric_center[0] = np.mean(x)
    geometric_center[1] = np.mean(y)
    geometric_center[2] = np.mean(z)

    # TODO: I am sure there should be tools that can do it directly for a mesh, considering the mesh panel size. We may look for that later

    return geometric_center



def get_points_from_gmsh_file(gmsh_file):
    '''
    Get points from GMSH file.

    Parameters
    ==========
    gmsh_file       : str, path of the gmsh file

    Return
    ======
    mesh_points     : ndarray,  points of the mesh. The last axis is for coordinates
    '''
    mesh = meshio.read(gmsh_file)
    mesh_points = mesh.points
    return mesh_points, mesh


def get_r_as_phi_theta_for_surface_mesh_points(mesh_points, theta_n_points, phi_n_points):
    '''
    # TODO: For a given array of points, and a meshgrid of theta and phi, we need to get r from the center of mass, so that r starts at the origin and ends at the mesh surface
    '''
    print('BUILD ME!')


def get_spherical_harmonics_for_mesh(gmsh_file, theta_n_points, phi_n_points):
    '''
    # TODO: Calculate spherical harmonics for a given mesh file
    '''
    # read points from the mesh file, OK!    
    mesh_points = get_points_from_gmsh_file(gmsh_file)

    # get r as a 2D array, a meshgrid in theta and phi
    # TODO
    r_2D_grid = get_r_as_phi_theta_for_surface_mesh_points(mesh_points, theta_n_points, phi_n_points)

    # calculate the coefficients based on geometry
    geom_based = pyshtools.SHGrid.from_array(np.array(r_2D_grid))
    #geom_based.lmax=5
    #geom_based.grid='DH'
    coeff_array = geom_based.expand(normalization='4pi', csphase=1, lmax_calc=12)

    print('BUILD ME!')



# Triangle mesh 














def tri_mesh(index,Tir_array,mesh_points):
    """
    contract tri_array and mesh_corr
    """
    Tir_array_test=(Tir_array[index])
    mesh_points_test=np.array([mesh_points[Tir_array_test[0]],mesh_points[Tir_array_test[1]],mesh_points[Tir_array_test[2]]])
    return Tir_array_test,mesh_points_test



    # Building triangle identifier
def tri_normal(tri_array,mesh_corr):
    '''
tri_array: input with 1D triangle array with vertex [v1,v2,v3]
mesh_corr: input with 3D arrary with the form [[],[],[]]. 
It represents the coresponding coordinate.

    '''
    # access vertex
    
    v1=tri_array[0]
    v2=tri_array[1]
    v3=tri_array[2]
    # Access datapoints
    v1_cor=mesh_corr[0]
    v2_cor=mesh_corr[1]
    v3_cor=mesh_corr[2]
    # Finding the area using cross product
        # Vectors
    V_i=v2_cor-v1_cor
    V_j=v3_cor-v1_cor

    Vector_area=np.cross(V_i,V_j)
        # Area
    Area=np.sqrt(Vector_area.dot(Vector_area))
        # unit direction
    dir_id=Vector_area/Area

        # geo center 
    center_tri=v1_cor#np.mean(mesh_corr,axis=0)

    return np.array([Area,dir_id,center_tri])
def dir_info(theta,phi):
    e_r=np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
    # [ [x...] ,[y...],[z....] ]
    return e_r

def tri_interp(Tri_info,dir_info):
    '''
    Tri_info gives the info from tri_normal
    dir_info gives the direction in an 1D array [theta,phi]

    '''
    # Normal vector for the plane 
    n_norm=Tri_info[:,0]*Tri_info[:,1]
    #print(n_norm)
    # vector const for plane 
    v_p0=Tri_info[:,2]
    # Unit vector dir 
    e_r=dir_info

    # Vector dot with const ( n dot v0=n dot v_desire)
    product_const=np.array(list(n_norm*v_p0))
    #print(product_const)
    # reshape
    n_norm=np.array(list(n_norm))
    product_coff=np.array(list(n_norm*e_r.T))

    const_dot=np.sum(product_const,axis=1)
    #print(const_dot)
    coff_dot=np.sum(product_coff,axis=1)
    #print(coff_dot)
    r_theta_phi=const_dot/coff_dot

    return r_theta_phi

def norm_id(tri_points,mesh_points):
    '''
    calculate norm for all triangles 
    '''
    # num of triangles
    loop_length=len(tri_points)
    Tig_norm=[]
    # array append 
    for i in range(loop_length):
        s=tri_mesh(i,tri_points,mesh_points)
        tri_array=s[0]
        mesh_corr=s[1]
        #Testing tri
        Tri_info=tri_normal(tri_array,mesh_corr)
        # correction for the outward direction 
        norm_dir=Tri_info[1]
        geo_center=Tri_info[2]
        if np.dot(geo_center,norm_dir)<0:
            Tri_info[1]=-1*Tri_info[1]
        Tig_norm.append(Tri_info)

    return np.array(Tig_norm)

def align_dir_tirangle(norm_array,target_dir):
    '''
    Find the tirangle that close to the target direction 
    norm_array: Info-array from norm_id
    target_dir: given theta and phi
    '''
    # transform to an location array for position vector 
    # Normalize
    Nor_dir_p_val=norm_array[:,2] 
    Dir_array_p=np.array(list(Nor_dir_p_val))
    # magnitude
    Dir_array_p_mag=np.sqrt(np.sum(Dir_array_p**2,axis=1))
    # boradcast to 3D array
    Dir_array_p_mag=np.stack((Dir_array_p_mag,Dir_array_p_mag,Dir_array_p_mag), axis=-1)
    Dir_array_p=Dir_array_p/(Dir_array_p_mag)

    # transform to an location array for position vector 
    Nor_dir_val=norm_array[:,1] 
    Dir_array=np.array(list(Nor_dir_val))
    # magnitude
    Dir_array_mag=np.sqrt(np.sum(Dir_array**2,axis=1))
    # boradcast to 3D array
    Dir_array_mag=np.stack((Dir_array_mag,Dir_array_mag,Dir_array_mag), axis=-1)
    Dir_array=Dir_array/(Dir_array_mag)


    
    # dot product respect to location

    dotproduct_dir_p=np.dot(Dir_array_p,target_dir)

    # dot product respect to normal direction
    dotproduct_dir=np.dot(Dir_array,target_dir)
    # total contri
    dot_total=dotproduct_dir_p+dotproduct_dir
    # find min value and index which is the most aligned
    index_dir=np.argmax(dot_total,axis=0)
    print(index_dir)
    #print(dotproduct_dir)
    #print(index_dir)
    value_dir=np.max(dotproduct_dir)
    #print(value_dir)
    #if np.size(index_dir)>1:
    #    index_dir=index_dir[0]
    #    value_dir=value_dir[0]
    #    print('critical angle, take the first value')

    return norm_array[index_dir]




        
        






