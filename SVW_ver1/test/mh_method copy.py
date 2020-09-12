# %% Importing
#########################################################
#########################################################

import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
# get path of this script
path_this_script = os.getcwd()
path_this_script = os.path.realpath(__file__)

# add the ./src/ path to the search path
path_this_script_splitted = os.path.split(path_this_script)
this_script_filename = path_this_script_splitted[1]
path_this_script_splitted = os.path.split(path_this_script_splitted[0])
path_to_src = os.path.join(path_this_script_splitted[0], 'src')
sys.path.append(path_to_src)  # I could have used sys.path.append('../src/'), but it didn't work with the debugger
path_to_cache = os.path.join(path_this_script_splitted[0], 'cache')
from SphTools import *



# %%
A=pd.read_csv('/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/simulation.csv')
B=A['simulation_result_filepath']

# %%
root_path='/home/melzouka/elzouka_codes_library/Python/SphEM/'
csv_path='/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/'
processed_data_dir='/home/zhw/SphEM-1/Processed_data'
path_to_cache = os.path.join(root_path, B[8])
get_data_from_SCUFF_TMATRIX_file(path_to_cache)


# %%

def pd_resize(pd_array,id_key,search_id='id'):
    pd_properties=pd_array.columns
   # print(pd_properties)
    arv_array=[]
    mis_array=[]
    for i in range(0,len(id_key)):
        check_id=pd_array[search_id]==id_key[i]
        s=pd_array[check_id]
        if s.empty:
            mis_array.append(id_key[i])
            
        arv_array.append(s)
       # print(s)
    
    arv_array=pd.concat(arv_array)
            
    return arv_array ,mis_array
    #print('--')
            #print(pd_array[pd_properties[i]])

def prop_finder(simu_type,prop,file_path):
    try:
        new_prop_array=[]
        if simu_type=='neq':
            NEQ_df=get_data_from_SCUFF_NEQ_file(file_path)

            try:
                prop_neq=NEQ_df.columns
                NEQ_df_prop=[]
                for j in range(0,len(prop)):
                    check_in=prop[j] in prop_neq 
                    #print(check_in)
                    if check_in:
                        NEQ_df_prop.append(NEQ_df[prop[j]])
                    else:
                        print('s')
                return NEQ_df_prop
            except: 
                pass

        elif simu_type=='scatter':
            Scatter_df=get_data_from_SCUFF_SCATTER_file(file_path)
            prop_scat=Scatter_df.columns
            try:
                Scat_df_prop=[]

                for j in range(0,len(prop)):
                    check_in=prop[j] in prop_scat 
                   # print(len(prop))
                   # print(check_in)
                    if check_in:
                        Scat_df_prop=Scatter_df[prop]
                    
                    else:
                        print('pp')
                return Scat_df_prop
            except:
                pass

        elif simu_type=='tmatrix':
            T_matrix_df=get_data_from_SCUFF_TMATRIX_file(file_path)
            
            prop_tmat=T_matrix_df.columns
            try:
                tmat_df_prop=[]
                for j in range(0,len(prop)):
                    check_in=prop[j] in prop_tmat 
                    #print(check_in)
                    if check_in:
                        tmat_df_prop=T_matrix_df[prop]
                        
                    else:
                        pass
                return tmat_df_prop
            except:
                pass

    except:
        raise TypeError

def path_gen(root_path,data_pd):
    file_path_ralative=data_pd['simulation_result_filepath'].to_numpy()
    sim_id=data_pd['simulation_id'].to_numpy()
    mesh_id=data_pd['mesh_id'].to_numpy()
    real_path_array=[]
    simu_id_array=[]
    mesh_id_array=[]
    for i in range(0,len(file_path_ralative)):
        #print(file_path_ralative[i])
        current_path=os.path.join(root_path, file_path_ralative[i])
        real_path_array.append(current_path)

    return real_path_array,sim_id,mesh_id












    
def search_id(mesh,Simulation,
                geom_id,mesh_id,material_id,prop_id):
    geo_id_array_at_mesh_file=pd_resize(pd_array=mesh,id_key=geom_id,search_id='geometry_id')[0]
    # check mesh id 
    mesh_id_array_at_mesh_file=pd_resize(pd_array=geo_id_array_at_mesh_file,id_key=mesh_id,search_id='mesh_id')
    # check simulation id 
    vali_mesh_id=mesh_id_array_at_mesh_file[0]['mesh_id'].to_numpy()

    simu_id_array_at_simu_file=pd_resize(pd_array=Simulation,id_key=vali_mesh_id,search_id='mesh_id')

    vali_simu_id=np.unique(simu_id_array_at_simu_file[0]['material_id'].to_numpy())
    material_id_array_at_simu_file=pd_resize(pd_array=simu_id_array_at_simu_file[0],id_key=vali_simu_id,search_id='material_id')
    print(material_id_array_at_simu_file)


    #print(vali_mesh_id)
    #print(vali_simu_id)
    simu_array_path=material_id_array_at_simu_file[0]['simulation_result_filepath']
    neq_type_simu=material_id_array_at_simu_file[0][material_id_array_at_simu_file[0]['simulation_type']=='neq']
    scatter_type_simu=material_id_array_at_simu_file[0][material_id_array_at_simu_file[0]['simulation_type']=='scatter']
    tmatrix_type_simu=material_id_array_at_simu_file[0][material_id_array_at_simu_file[0]['simulation_type']=='tmatrix']
    neq_p=path_gen(root_path,neq_type_simu)
    scatter_p=path_gen(root_path,scatter_type_simu)
    tmatrix_p=path_gen(root_path,tmatrix_type_simu)
    target_path_neq=neq_p[0]
    target_path_scatter=scatter_p[0]
    target_path_tmatrix=tmatrix_p[0]

    scatter_mesh_id=scatter_p[2]
    scatter_simu_id=scatter_p[1]
   

    #print(len(target_path))
    #print(neq_type_simu['simulation_result_filepath'])
    #print(neq_type_simu)
    neq_prop=[]
    for n in range(0,len(target_path_neq)):
        file_path=target_path_neq[n]
        K=prop_finder('neq',prop_id,file_path)
        neq_prop.append(K)
        
    scatter_prop=[]
    for n in range(0,len(target_path_scatter)):
        file_path=target_path_scatter[n]
        
        K=prop_finder('scatter',prop_id,file_path)

        scatter_prop.append(K)
    #scatter_prop=pd.concat(scatter_prop,axis=1)
    #id_ver_scatter_prop=pd.concat([scatter_prop,scatter_mesh_id,scatter_simu_id],axis=1)
    


    tmatrix_prop=[]
    for n in range(0,len(target_path_tmatrix)):
        file_path=target_path_tmatrix[n]

        K=prop_finder('tmatrix',prop_id,file_path)
        tmatrix_prop.append(K)
    print(scatter_prop)

    
    #simu_type=simu_type
    #for i in range(0,len(vali_mesh_id)):
    #    prop_finder(simu_type,prop_name=0,file_path)
    
    #except:
        
    return mesh
#%%
geom_id=['wireSphEnds/Lx_0.3/D_0.04/','wire/Lx_0.4/D_0.08/','spheroid/rx_0.04/rz_0.08/']
mesh_id=['wireSphEnds/Lx_0.3/D_0.04//clscale_1_clmax_1','wireSphEnds/Lx_0.3/D_0.04//clscale_0.95_clmax_1.2','spheroid2/rx_0.04/rz_0.08//clscale_0.95_clmax_1','adhihadih']
simu_id=[2,3,4,5,34,1,2,4,5]
material_id=['CONST_EPS_1.8//Au']
root_path='/home/melzouka/elzouka_codes_library/Python/SphEM/'
csv_path='/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/'
processed_data_dir='/home/zhw/SphEM-1/Processed_data'
prop_id=['scatter_CrossSection_m2','freq_radpers']
mesh=F[0]
Simulation=F[1]
AA=search_id(mesh,Simulation,geom_id,mesh_id,material_id,prop_id)
file_look(mesh,Simulation,geom_id)
#%%

class pointer():
    def __init__(self,root_path,csv_path,
                    processed_data_dir):
        self.root_path=root_path
        self.csv_path=csv_path
        self.processed_data_dir=processed_data_dir
        sim_csv_path=os.path.join(csv_path,  'simulation.csv')
        geo_csv_path=os.path.join(csv_path, 'geometry.csv')
        mater_csv_path=os.path.join(csv_path, 'material.csv')
        msh_csv_path=os.path.join(csv_path, 'mesh.csv')

        Simulation=pd.read_csv(sim_csv_path)
        self.Geometry=pd.read_csv(geo_csv_path)
        self.materials=pd.read_csv(mater_csv_path)
        self.mesh=pd.read_csv(msh_csv_path)
        Simulation["material_id"] = Simulation["material_surr"] +'//'+ Simulation["material_particle"]
        self.Simulation=Simulation

        self.geo_id_list=np.array(self.mesh['geometry_id'])
        self.mesh_id_list=np.array(self.mesh['mesh_id'])
        self.simu_id_list=np.array(self.Simulation['simulation_id'])
        self.material_id_list=np.array(self.Simulation['material_id'])
        print(self.Simulation.columns)

    def meshsimu(self):    
        return self.mesh, self.Simulation

    def get_id(self):
        id_list={'geo_id_list':self.geo_id_list,
        'mesh_id_list':self.mesh_id_list, 
        'simu_id_list':self.simu_id_list,'material_id_list':self.material_id_list}
        import csv
        with open(os.path.join(self.processed_data_dir,  'id_list.csv'), 'w') as f:  # Just use 'w' mode in 3.x
            w = csv.DictWriter(f, id_list.keys())
            w.writeheader()
            w.writerow(id_list)
        return id_list

  
        #print(geo_id_list)
    def data_init(self, geom_id,mesh_id,simu_type,prop_id):
        # getting all files
        root_path=self.root_path
        csv_path=self.csv_path
        processed_data_dir=self.processed_data_dir

        Simulation=self.Simulation
        Geometry=self.Geometry
        materials=self.materials
        mesh=self.mesh

        geo_check=len(geom_id)
        mesh_check=len(mesh_id)
        simu_check=len(simu_id)
        if  (geo_check<=mesh_check):
            # check geo id
            #try:
            geo_id_array_at_mesh_file=pd_resize(pd_array=mesh,id_key=geom_id,search_id='geometry_id')[0]
            # check mesh id 
            mesh_id_array_at_mesh_file=pd_resize(pd_array=geo_id_array_at_mesh_file,id_key=mesh_id,search_id='mesh_id')
            # check simulation id 
            vali_mesh_id=mesh_id_array_at_mesh_file[0]['mesh_id'].to_numpy()

            simu_id_array_at_simu_file=pd_resize(pd_array=Simulation,id_key=vali_mesh_id,search_id='mesh_id')

            #print(simu_id_array_at_simu_file)
            simu_array_path=simu_id_array_at_simu_file[0]['simulation_result_filepath']
            neq_type_simu=simu_id_array_at_simu_file[0][simu_id_array_at_simu_file[0]['simulation_type']=='neq']
            scatter_type_simu=simu_id_array_at_simu_file[0][simu_id_array_at_simu_file[0]['simulation_type']=='scatter']
            tmatrix_type_simu=simu_id_array_at_simu_file[0][simu_id_array_at_simu_file[0]['simulation_type']=='tmatrix']
            neq_p=path_gen(self.root_path,neq_type_simu)
            scatter_p=path_gen(self.root_path,scatter_type_simu)
            tmatrix_p=path_gen(self.root_path,tmatrix_type_simu)
            target_path_neq=neq_p[0]
            target_path_scatter=scatter_p[0]
            target_path_tmatrix=tmatrix_p[0]
            print(target_path_tmatrix)

            scatter_mesh_id=scatter_p[2]
            scatter_simu_id=scatter_p[1]
            

            #print(len(target_path))
            #print(neq_type_simu['simulation_result_filepath'])
            #print(neq_type_simu)
            neq_prop=[]
            for n in range(0,len(target_path_neq)):
                file_path=target_path_neq[n]
                K=prop_finder('neq',prop_id,file_path)
                neq_prop.append(K)
                
            scatter_prop=[]
            for n in range(0,len(target_path_scatter)):
                file_path=target_path_scatter[n]
                
                K=prop_finder('scatter',prop_id,file_path)
 
                scatter_prop.append(K)
            #scatter_prop=pd.concat(scatter_prop,axis=1)
            #id_ver_scatter_prop=pd.concat([scatter_prop,scatter_mesh_id,scatter_simu_id],axis=1)
            


            tmatrix_prop=[]
            for n in range(0,len(target_path_tmatrix)):
                file_path=target_path_tmatrix[n]
  
                K=prop_finder('tmatrix',prop_id,file_path)
                tmatrix_prop.append(K)
            #print(scatter_prop)
    
            
            #simu_type=simu_type
            #for i in range(0,len(vali_mesh_id)):
            #    prop_finder(simu_type,prop_name=0,file_path)
            
            #except:
                
            return mesh

        elif geo_check==1 and (not mesh_check) and ( not simu_check):
            pass
        elif geo_check==1 and (geo_check<=mesh_check) and ( not simu_check):
            pass
        elif (not geo_check) and (not mesh_check) and ( simu_check):
            pass
        elif (not geo_check) and (mesh_check) and ( not simu_check):
            pass


        
        
        #path_to_cache = os.path.join(root_path, B[8])
     

    #def add_branch(self,)



# %%
R=pointer(root_path,csv_path,processed_data_dir)
geom_id=['wireSphEnds/Lx_0.3/D_0.04/','wire/Lx_0.4/D_0.08/','spheroid/rx_0.04/rz_0.08/']
mesh_id=['wireSphEnds/Lx_0.3/D_0.04//clscale_1_clmax_1','wireSphEnds/Lx_0.3/D_0.04//clscale_0.95_clmax_1.2','spheroid2/rx_0.04/rz_0.08//clscale_0.95_clmax_1','adhihadih']
simu_id=[2,3,4,5,34,1,2,4,5]
material_id=['CONST_EPS_1.8//Au']
root_path='/home/melzouka/elzouka_codes_library/Python/SphEM/'
csv_path='/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/'
processed_data_dir='/home/zhw/SphEM-1/Processed_data'
prop_id=['scatter_CrossSection_m2','freq_radpers']
K=R.data_init(geom_id,mesh_id,simu_id,prop_id)
F=R.meshsimu()

# %%
file_csv_path='/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/geometry.csv'
master_name='geometry'

# %%
R.get_id()

# %%


# %%
# a function identify files
def pd_resize_n(pd_array,id_key,search_id='id'):
    pd_properties=pd_array.columns
   # print(pd_properties)
    arv_array=[]
    mis_array=[]
    for i in range(0,len(id_key)):
        check_id=pd_array[search_id]==id_key[i]
        s=pd_array[check_id]
        if s.empty:
            mis_array.append(id_key[i])
            
        arv_array.append(s)
       # print(s)
    
    arv_array=pd.concat(arv_array)
    '''
    neq_type_simu=arv_array[arv_array['simulation_type']=='neq']
    scatter_type_simu=arv_array[arv_array['simulation_type']=='scatter']
    tmatrix_type_simu=arv_array[arv_array['simulation_type']=='tmatrix']

    neq_full_path=path_gen(root_path,neq_type_simu)
    scatter_full_path=path_gen(root_path,scatter_type_simu)
    tmatrix_full_path=path_gen(root_path,tmatrix_type_simu)
    print(reduced_files)
    '''
    return arv_array ,mis_array
def full_path(root_path,relative_path_contained_data_array):
    check_len=len(relative_path_contained_data_array)
    
    relative_path_contained_data_array['full_path']=root_path+relative_path_contained_data_array['simulation_result_filepath']
    return relative_path_contained_data_array

def full_path_concat(full_path_data_array):
    # find path 
    neq_type_simu=full_path_data_array[full_path_data_array['simulation_type']=='neq']
    scatter_type_simu=full_path_data_array[full_path_data_array['simulation_type']=='scatter']
    tmatrix_type_simu=full_path_data_array[full_path_data_array['simulation_type']=='tmatrix']
     
    return full_path_data_array,neq_type_simu, scatter_type_simu,tmatrix_type_simu


    #except:
      #  print('s')
def file_look(root_path,mesh,Simulation,geom_id=[],mesh_id=[],material_id=[],simulation_id=[]):
    # getting the size of the input
    geo_size=len( geom_id)
    mesh_szie=len(mesh_id)
    material_size=len(material_id)
    simulation_size=len(simulation_id)

    # check from back to start 

    if simulation_size==0:
        No_found_geo_id=[]
        No_found_simu_id=[]
        No_found_material_id=[]
        No_found_mesh_id=[]
        # check material and mesh size
        if  mesh_szie!=0:
            # look mesh id
            Full_files_log=pd_resize_n(Simulation,mesh_id,search_id='mesh_id')

            reduced_files=Full_files_log[0]
            No_found_mesh_id=Full_files_log[1]
            No_found_geo_id=[]
            No_found_simu_id=[]
            No_found_material_id=[]
            # path gene
            full_path_array=full_path(root_path,reduced_files)
            # file display
            



            if material_size!=0:
                # checking material from the reduced array
                Full_files_log=pd_resize_n(full_path_array,material_id,search_id='material_id')
                #print(Full_files_log)
                material_reduced_files=Full_files_log[0]
                No_found_material_id=Full_files_log[1]

                Selected_property_files=full_path_concat(material_reduced_files)# (NEQ, Scatter, Tmatrix)

            else:
                Selected_property_files=full_path_concat(full_path_array)
                No_found_material_id=[]

            # print out the log file
            full_path_data_array=Selected_property_files[0]
            return_file_log=pd.DataFrame([
                full_path_data_array['simulation_id'],
                full_path_data_array['material_id'],
                full_path_data_array['simulation_type'],
                full_path_data_array['full_path']
            ]).T
             

        elif mesh_szie==0 and material_size!=0:
            Full_files_log=pd_resize_n(Simulation,material_id,search_id='material_id')
            reduced_files=Full_files_log[0]
            No_found_material_id=Full_files_log[1]
            No_found_mesh_id=[]

            # path gene
            full_path_array=full_path(root_path,reduced_files)  
            Selected_property_files=full_path_concat(full_path_array) 

            # print out the log file
            full_path_data_array=Selected_property_files[0]
            return_file_log=pd.DataFrame([
                full_path_data_array['simulation_id'],
                full_path_data_array['material_id'],
                full_path_data_array['simulation_type'],
                full_path_data_array['full_path']
            ]).T  
            

        elif mesh_szie==0 and material_size==0 and geo_size!=0:
            # checking file
            Full_files_log_mesh=pd_resize_n(mesh,geom_id,search_id='geometry_id')
            reduced_files_mesh_from_geo=Full_files_log_mesh[0]
            No_found_mesh_id=[]
            No_found_material_id=[]
            No_found_geo_id=Full_files_log_mesh[1]

             
            # grab mesh id from geo
            mesh_id_from_geo=reduced_files_mesh_from_geo['mesh_id'].to_numpy()

            # Checking new mesh file. Same process from above
            Full_files_log=pd_resize_n(Simulation,mesh_id_from_geo,search_id='mesh_id')
            reduced_files=Full_files_log[0]
            No_found_mesh_id=Full_files_log[1]
            # path gene
            full_path_array=full_path(root_path,reduced_files)
            Selected_property_files=full_path_concat(full_path_array) 

            # print out the log file
            full_path_data_array=Selected_property_files[0]
            return_file_log=pd.DataFrame([
                full_path_data_array['simulation_id'],
                full_path_data_array['material_id'],
                full_path_data_array['simulation_type'],
                full_path_data_array['full_path']
            ]).T  

        else:
            Selected_property_files=([],[],[],[])
            No_found_mesh_id=[]
            No_found_material_id=[]
            No_found_geo_id=[]
            

    else:
            Full_files_log=pd_resize_n(Simulation,simulation_id,search_id='simulation_id')
            reduced_files=Full_files_log[0]
            No_found_simu_id=Full_files_log[1]
            # path gene
            full_path_array=full_path(root_path,reduced_files)
            Selected_property_files=full_path_concat(full_path_array) 

            # print out the log file
            full_path_data_array=Selected_property_files[0]
            return_file_log=pd.DataFrame([
                full_path_data_array['simulation_id'],
                full_path_data_array['material_id'],
                full_path_data_array['simulation_type'],
                full_path_data_array['full_path']
            ]).T          
            No_found_mesh_id=[]
            No_found_geo_id=[]
            No_found_material_id=[]

            #print(Selected_property_files)
 
        

    return Selected_property_files,return_file_log


# Checking the properties of three files 
def full_properties_look(neq_pandas,
                    scatter_pandas,
                    tmatrix_pandas):
    # get the file path 
    
    neq_file_path=neq_pandas['full_path'].reset_index(drop=True)
    neq_len=len(neq_file_path)
    #neq_prop_len=len(neq_properties)
    neq_simu_id=neq_pandas['simulation_id'].reset_index(drop=True)

    
    scatter_file_path=scatter_pandas['full_path'].reset_index(drop=True)
    scatter_len=len(scatter_file_path)
    #scatter_prop_len=len(scatter_properties)
    scatter_simu_id=scatter_pandas['simulation_id'].reset_index(drop=True)

    
    tmatrix_file_path=tmatrix_pandas['full_path'].reset_index(drop=True)
    tmatrix_len=len(tmatrix_file_path)
    #tmatrix_prop_len=len(tmatrix_properties)
    tmatrix_simu_id=tmatrix_pandas['simulation_id'].reset_index(drop=True)
 
    # EXPAND the file
    if neq_len!=0:
        NEQ_properties=[]
        for i in range(0,neq_len):
            # define the current simu_id
            current_simu_path=neq_file_path[i]
            current_simu_id=neq_simu_id[i]
            Current_NEQ_file=get_data_from_SCUFF_NEQ_file(current_simu_path)
            Current_NEQ_file['simulation_id']=current_simu_id
            NEQ_properties.append(Current_NEQ_file)
            # combine into one pandas frame 
        NEQ_properties=pd.concat(NEQ_properties)
    else:
        NEQ_properties=pd.DataFrame()


    if scatter_len!=0:
        Scatter_properties=[]
        for i in range(0,scatter_len):
            
            # define the current simu_id
            current_simu_path=scatter_file_path[i]
            current_simu_id=scatter_simu_id[i]
            Current_SCATTER_file=get_data_from_SCUFF_SCATTER_file(current_simu_path)
            Current_SCATTER_file['simulation_id']=current_simu_id
            Scatter_properties.append(Current_SCATTER_file)
            # combine into one pandas frame 
        Scatter_properties=pd.concat(Scatter_properties)

    else:
        Scatter_properties=pd.DataFrame()

    if tmatrix_len!=0:
        Tmatrix_properties=[]
        for i in range(0,scatter_len):
            
            # define the current simu_id
            current_simu_path=tmatrix_file_path[i]
            current_simu_id=tmatrix_simu_id[i]
            Current_TMATRIX_file=get_data_from_SCUFF_TMATRIX_file(current_simu_path)
            Current_TMATRIX_file['simulation_id']=current_simu_id
            Tmatrix_properties.append(Current_TMATRIX_file)
            # combine into one pandas frame 
        Tmatrix_properties=pd.concat(Tmatrix_properties)

        #print(Tmatrix_properties)
    else:
        Tmatrix_properties=pd.DataFrame()

    return NEQ_properties, Scatter_properties, Tmatrix_properties
 

def properties_look(NEQ_properties, neq_properties,
                    Scatter_properties,scatter_properties,
                    Tmatrix_properties,tmatrix_properties):
    # check is empty 
    neq_len=len(NEQ_properties)
    scatter_len=len(Scatter_properties)
    tmatrix_len=len(Tmatrix_properties)

    try:
        if neq_len!=0:
            # append simu id 
            neq_properties+=['simulation_id']
            Target_NEQ=NEQ_properties[neq_properties]
            Target_NEQ.set_index(['simulation_id'],inplace=True)
        else:
            Target_NEQ=pd.DataFrame()
    except:
        Target_NEQ=['wrong properties name']

    try:
        if scatter_len!=0:
            # append simu id 
            
            scatter_properties+=['simulation_id']
    
            Target_Scatter=Scatter_properties[scatter_properties]
            Target_Scatter.set_index(['simulation_id'],inplace=True)
        else:
            Target_Scatter=pd.DataFrame()
    except:
        Target_Scatter=['wrong properties name']

    try:
        if tmatrix_len!=0:
            # append simu id 
            tmatrix_properties+=['simulation_id']
            Target_Tmatrix=Tmatrix_properties[tmatrix_properties]
            Target_Tmatrix.set_index(['simulation_id'],inplace=True)
        else:
            Target_Tmatrix=pd.DataFrame()
    except:
        Target_Tmatrix=['wrong properties name']

    return Target_NEQ, Target_Scatter, Target_Tmatrix

geom_id=[]#['wireSphEnds/Lx_0.3/D_0.04/','wire/Lx_0.4/D_0.08/','spheroid/rx_0.04/rz_0.08/']
mesh_id=['wireSphEnds/Lx_0.3/D_0.04//clscale_1_clmax_1','wireSphEnds/Lx_0.3/D_0.04//clscale_0.95_clmax_1.2','spheroid/rx_0.04/rz_0.08//clscale_0.95_clmax_1','adhihadih']
simulation_id=[]#['wireSphEnds/Lx_0.3/D_0.04//clscale_0.9_clmax_1.2/CONST_EPS_1.8/Au/scatter/tol_1e-05']#[2,3,4,5,34,1,2,4,5]
material_id=['CONST_EPS_1.8//Au']
mesh=F[0]
Simulation=F[1]

S=file_look(root_path,mesh,Simulation,geom_id,mesh_id,material_id,simulation_id)
Prop=S[0]
neq_pandas=Prop[1]
scatter_pandas=Prop[2]
tmatrix_pandas=Prop[3]
neq_prop=[]
scatter_prop=['scatter_CrossSection_m2','freq_radpers']
tmatrix_prop=[]
P=full_properties_look(neq_pandas,
                    scatter_pandas,
                    tmatrix_pandas)
NEQ_properties=P[0]
Scatter_properties=P[1]
Tmatrix_properties=P[2]

SS=properties_look(NEQ_properties, neq_prop,
                    Scatter_properties,scatter_prop,
                    Tmatrix_properties,tmatrix_prop)
# %%




# %%
