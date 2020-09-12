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


#%%

def alpha2array(Chart):
    Chart=''.join(filter(str.isalpha, Chart))
    Alpha_basis=np.eye(26)
    
    try:
        chart=list(Chart.lower())
        
        
        for i in range(0,len(chart)):
            current_alpha=chart[i]

            current_index_alpha=ord(current_alpha)-97
            current_array_alpha=Alpha_basis[current_index_alpha,:]
            if i==0:
                alpha_array=current_array_alpha
            else:
                alpha_array=alpha_array+current_array_alpha
            
        
        str_array=(''.join(str(int(x))for x in (alpha_array)))
 
        char={Chart:str_array }
        return alpha_array,char

    except:
        alpha_array=np.ones((26,))
        if Chart=='':
            alpha_array=np.zeros((26,))
        str_array=(''.join(str(int(x))for x in (alpha_array)))
        char={str(Chart):str_array }
        return alpha_array,char

def float2array(Float):
    Alpha_basis=np.eye(26)
    # init 
    try: 
        float_init=float(str(Float))
        for i in range(0,26):
            current_float_init=int(str(float_init/26)[-1])

            if i==0:
                alpha_array=Alpha_basis[i,:]*current_float_init

            else:
                alpha_array=alpha_array+Alpha_basis[i,:]*current_float_init

            float_init=float_init/26

        str_array=(''.join(str(int(x))for x in (alpha_array)))

        flot={str(Float):str_array }
        return alpha_array,flot

    except:
        if str(Float)=='nan':
            alpha_array=np.zeros((26,))

            str_array=(''.join(str(int(x))for x in (alpha_array)))
            flot={str(Float):str_array }
            return alpha_array,flot
        else:
            raise TypeError
            
def Merge_dic(dict1, dict2): 
    return(dict2.update(dict1))


class pointing_device():
    import os 
    import sys
    import pandas as pd
    import numpy as np

        
    def __init__(self):
        pass 




    def init_data_base(self,root_path,
    csv_path,processed_data_dir):
        A=pd.read_csv(csv_path)
        
        path_to_cache = os.path.join(root_path, B[8])
        self.raw_data_dir= raw_data_dir
        raw_file_list=[]
        for root, dirs, files in os.walk(raw_data_dir):
        
            print(files[0])
        
        
        #print(files)
        Current_file_full_path=os.path.join(raw_data_dir, files[0])
        S=pd.read_pickle(Current_file_full_path)
        Atr=S.columns
        print(S.columns)
        current_set_index=1
        
        for i in range(1,len(Atr)):
            Current_index=i
            atr=Atr[Current_index]
            #Attribute
            Current_art_id=alpha2array(atr)
            Current_atr_data=S[atr]
            data_type=Current_atr_data.dtype
            print(atr)
            Current_atr_data=np.array(Current_atr_data)
            try:
                if (data_type=='int64') or (data_type=='float64'):
                    Current_id=float2array(Current_atr_data[current_set_index])

                if data_type=='object':
                    check_str=Current_atr_data[current_set_index]
                    #print(check_str)
                    if isinstance(check_str, str)==True:
                        Current_id=alpha2array(check_str)
                    else: 
                        Current_id=(np.ones((26,))*1)          
                        str_array=(''.join(str(int(x))for x in (Current_id)))
                        Current_id_str={str(check_str):str_array }
                        Current_id=(Current_id,Current_id_str)

                        # expand the nested array 
                        Current_data=Current_atr_data[current_set_index]
                        
                        try: 
                            Data_columns=Current_data.columns
                            print(Data_columns)
                            for i in range(0,len(Data_columns)):
                                Current_data_field=Current_data[Data_columns[i]]
                                Current_data_type_field=Current_data_field.dtype
                                print(Current_data_type_field)
                               # try:
                               #     if (Current_data_type_field=='int64') or (Current_data_type_field=='float64'):
                                       # Current_id_field=float2array(Current_atr_data[current_set_index])

                              #      if Current_data_type_field=='object':
                                    #    check_str=Current_atr_data[current_set_index]




                        except:
                            if isinstance(np.sum(Current_data),float):
                                print('yes')
                            else:
                                pass 
                        #print((Current_data))
                            #print(Current_data)
                         



            except:
                Current_id=(np.ones((26,))*-1)          
                str_array=(''.join(str(int(x))for x in (Current_id)))
                Current_id_str={str(check_str):str_array }
                Current_id=(Current_id,Current_id_str)

                print('error_naming')
            #print(Current_id)
            if i==1:
                #print(Current_id)
                
                ID_matrix=Current_id[0]
            else:
                #print(atr)
                #print(Current_id[0].shape)
                ID_matrix=np.vstack([ID_matrix,Current_id[0]])
                #print(Current_id)

        #print(ID_matrix)


        return S



            #print(Current_atr_data)



         #= os.path.join(raw_data_dir, 'src')
    #    all_data_orig=pd.read_pickle('/home/melzouka/elzouka_codes_library/Python/SphEM/results/digested_for_Steve/all_data_tmatrix_for_Steve_fixed_tmatrix.pkl')
   #     try:
  #          Defalt_data_dir=os.mkdir('/home/zhw/SphEM-1/Processed_data')
  #      except:
  #          os.rmdir('/home/zhw/SphEM-1/Processed_data')
      #      Defalt_data_dir=os.mkdir('/home/zhw/SphEM-1/Processed_data')
        
       # print(processed_data_dir)




#%%
import pandas as pd
import numpy as np
A=pointing_device( )
R=A.init_data_base(root_path='/home/melzouka/elzouka_codes_library/Python/SphEM/',csv_path='/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/simulation.csv',processed_data_dir='/home/zhw/SphEM-1/Processed_data')

# %%
s=alpha2array('')


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
class pointer():
    def __init__(self):
        pass
    def data_init(self,root_path,csv_path,processed_data_dir):
        sim_csv_path=os.path.join(csv_path,  'simulation.csv')
        geo_csv_path=os.path.join(csv_path, 'geometry.csv')
        mater_csv_path=os.path.join(csv_path, 'material.csv')
        msh_csv_path=os.path.join(csv_path, 'mesh.csv')

        Simulation=pd.read_csv(sim_csv_path)
        Geometry=pd.read_csv(geo_csv_path)
        materials=pd.read_csv(mater_csv_path)
        mesh=pd.read_csv(msh_csv_path)
        #print(mesh)
        #path_to_cache = os.path.join(root_path, B[8])
     
    def add_master(self,file_csv_path,master_name):
        try:
            Master_data_array=pd.read_csv(file_csv_path)
            master=np.array(Master_data_array[master_name])
            
            for i in range(0,len(master)):
                current_master_id=alpha2array(master[i])
                if i==0:
                    dict1=current_master_id[1]
                    
                    
                    master_id=current_master_id[0]
                else:
                    dict1={**dict1 , **current_master_id[1]}
                    #print(current_master_id[1])
                    master_id=np.vstack([master_id,current_master_id[0]])
                    print(dict1)


            print(master)
        except:
            print('ee')

    #def add_branch(self,)



# %%
R=pointer()
R.data_init(root_path='/home/melzouka/elzouka_codes_library/Python/SphEM/',csv_path='/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/',processed_data_dir='/home/zhw/SphEM-1/Processed_data')

# %%
file_csv_path='/home/melzouka/elzouka_codes_library/Python/SphEM/results/database/scratch/geometry.csv'
master_name='geometry'
R.add_master(file_csv_path,master_name)

# %%
A-pd.