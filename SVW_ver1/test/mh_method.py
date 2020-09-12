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

def search_id(mesh,simu,geom_id,mesh_id,material_id=0,simu_type,prop_id):
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

        scatter_mesh_id=scatter_p[2]
        scatter_simu_id=scatter_p[1]
        print(scatter_mesh_id)

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

            scatter_mesh_id=scatter_p[2]
            scatter_simu_id=scatter_p[1]
            print(scatter_mesh_id)

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
root_path='/home/melzouka/elzouka_codes_library/Python/SphEM/'
csv_path='/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/'
processed_data_dir='/home/zhw/SphEM-1/Processed_data'
prop_id=['scatter_CrossSection_m2','freq_radpers']
K=R.data_init(geom_id,mesh_id,simu_id,prop_id)

# %%
file_csv_path='/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/geometry.csv'
master_name='geometry'

# %%
R.get_id()

# %%

 


# %%
