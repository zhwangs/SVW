# %% Importing
#########################################################
#########################################################

import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
# get path of this script
SphTools_full_path='/home/melzouka/elzouka_codes_library/Python/SphEM/src' # MAHMOUD: I have changed this path to make it pointing to the most updated file. # '/home/zhw/SphEM-1/src'
sys.path.append(SphTools_full_path)
from SphTools import *



# MAHMOUD: general notes:
'''
- Thanks Steve for making the class!

Suggestions for the code:
    - I see you have a method 'file_look' that can return simulation_ids given one or more of [geom_id,mesh_id,material_id,simu_id].
        I like it and it is very flexible.
        I think it would be helpful to add more functionality to your class, +by adding a method that return simulation_ids
        for the simulations that have exaclty the same parameters, except for one.
        -- for example, we have here the following column parameters that define each simulation:
            ['simulation_id', 'simulation_type', 'mesh_id', 'material_surr', 'material_particle', 'tolerance', 'RelOrient_id']
            I need to extract the simulation_ids for clusters of simuation that agrees on all columns, except for one column.

                example: I need to get simulation_ids for simulations that differs only in 'mesh_id', 
                so I need to group the DataFrame using all columns except for 'mesh_id'
                Then, I return to user only groups with 2 or more simulation_ids.

                You may use this function, or copy and modify it in any way you prefer:
                   SphTools > investigate_ind_dep_parameters  

            method(columnn_title_I_need_to_be_varied):

            columnn_title_I_need_to_be_varied: can be any one of the following:
            ['simulation_type', 'mesh_id', 'material_surr', 'material_particle', 'tolerance', 'RelOrient_id']

            ADVANCED:
            can we add columns from geometry.csv, mesh.csv





        -- I hope we have a visualization methods in our class, but I am not sure which is the best. Random ideas are here:
            --- given simulation_ids, plot a specific property (scattering_cross_section) VS another property (wavelength or angle) on the same plot, 
                and displaya ;egend that shows how the multiple lines in this plot differs (for example, they can have everything the same, except for the diameter 'D'. So, it would be helpful of we add a legend that contains only 'D = VALUE')

            --- given simulation_ids, plot (2D, 3D, or whatever convenient) how the simulations are agreeing or different. This would be very helpful for PIs to check on the progress of building our dataset.

        It would be really helpful if you may show examples on how to use each method.

LOW PRIORITY:
    - formatting to enhance code readability:
        -- it is preferred to create a space after each comma ', '
        -- space before and after the equal sign ' = '
    - 
'''



#%%

class pointer( ):
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


    def meshsimu(self):    
        return self.mesh, self.Simulation
# a function identify files
    def pd_resize_n(self,pd_array,id_key,search_id='id',check_init=False):
        arv_array=[]
        mis_array=[]
        if check_init==True:
            for i in range(0,len(id_key)):
                # check contain 
                check_id=pd_array[search_id].str.contains(id_key[i], regex=False)
                s=pd_array[check_id]
                if s.empty:
                    mis_array.append(id_key[i])
                    
                arv_array.append(s)
            
            arv_array=pd.concat(arv_array)                
        else:
            for i in range(0,len(id_key)):
                check_id=pd_array[search_id]==id_key[i]
                s=pd_array[check_id]
                if s.empty:
                    mis_array.append(id_key[i])
                    
                arv_array.append(s)
            
            arv_array=pd.concat(arv_array)
        return arv_array ,mis_array
# A function returns full path
    def full_path(self,root_path,relative_path_contained_data_array):
        check_len=len(relative_path_contained_data_array)
        
        relative_path_contained_data_array['full_path']=root_path+relative_path_contained_data_array['simulation_result_filepath']
        return relative_path_contained_data_array
# pandas with a full path column     
    def full_path_concat(self,full_path_data_array):
        # find path 
        neq_type_simu=full_path_data_array[full_path_data_array['simulation_type']=='neq']
        scatter_type_simu=full_path_data_array[full_path_data_array['simulation_type']=='scatter']
        tmatrix_type_simu=full_path_data_array[full_path_data_array['simulation_type']=='tmatrix']
        
        return full_path_data_array,neq_type_simu, scatter_type_simu,tmatrix_type_simu

# find target file 
    def file_look(self,geom_id=[],mesh_id=[],material_id=[],simulation_id=[]):
        root_path=self.root_path
        mesh=self.mesh
        Simulation=self.Simulation
        
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
                Full_files_log=self.pd_resize_n(Simulation,mesh_id,search_id='mesh_id')

                reduced_files=Full_files_log[0]
                No_found_mesh_id=Full_files_log[1]
                No_found_geo_id=[]
                No_found_simu_id=[]
                No_found_material_id=[]
                # path gene
                full_path_array=self.full_path(root_path,reduced_files)
                # file display
                



                if material_size!=0:
                    # checking material from the reduced array
                    Full_files_log=self.pd_resize_n(full_path_array,material_id,search_id='material_id')
                    material_reduced_files=Full_files_log[0]
                    No_found_material_id=Full_files_log[1]

                    Selected_property_files=self.full_path_concat(material_reduced_files)# (NEQ, Scatter, Tmatrix)
                else:
                    Selected_property_files=self.full_path_concat(full_path_array)
                    No_found_material_id=[]

                # print out the log file
                full_path_data_array=Selected_property_files[0]
                return_file_log=pd.DataFrame([
                    full_path_data_array['simulation_id'],
                    full_path_data_array['material_id'],
                    full_path_data_array['simulation_type'],
                    full_path_data_array['full_path']
                ]).T
                

            elif mesh_szie==0 and material_size!=0 and geo_size==0:
                Full_files_log=self.pd_resize_n(Simulation,material_id,search_id='material_id')
                reduced_files=Full_files_log[0]
                No_found_material_id=Full_files_log[1]
                No_found_mesh_id=[]

                # path gene
                full_path_array=self.full_path(root_path,reduced_files)  
                Selected_property_files=self.full_path_concat(full_path_array) 

                # print out the log file
                full_path_data_array=Selected_property_files[0]
                return_file_log=pd.DataFrame([
                    full_path_data_array['simulation_id'],
                    full_path_data_array['material_id'],
                    full_path_data_array['simulation_type'],
                    full_path_data_array['full_path']
                ]).T  
                
            elif mesh_szie==0 and material_size!=0 and geo_size!=0:
                Full_files_log=self.pd_resize_n(Simulation,material_id,search_id='material_id')
                reduced_files_material=Full_files_log[0]
                No_found_material_id=Full_files_log[1]
                No_found_mesh_id=[]
                # get mesh id, and check with geo_id 
                Full_files_log=self.pd_resize_n(reduced_files_material,geom_id,search_id='mesh_id',check_init=True)
                reduced_files=Full_files_log[0]
                No_found_geo_id=Full_files_log[1]

                
                # path gene
                full_path_array=self.full_path(root_path,reduced_files)  
                Selected_property_files=self.full_path_concat(full_path_array) 

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
                Full_files_log_mesh=self.pd_resize_n(mesh,geom_id,search_id='geometry_id')
                reduced_files_mesh_from_geo=Full_files_log_mesh[0]
                No_found_mesh_id=[]
                No_found_material_id=[]
                No_found_geo_id=Full_files_log_mesh[1]

                
                # grab mesh id from geo
                mesh_id_from_geo=reduced_files_mesh_from_geo['mesh_id'].to_numpy()

                # Checking new mesh file. Same process from above
                Full_files_log=self.pd_resize_n(Simulation,mesh_id_from_geo,search_id='mesh_id')
                reduced_files=Full_files_log[0]
                No_found_mesh_id=Full_files_log[1]
                # path gene
                full_path_array=self.full_path(root_path,reduced_files)
                Selected_property_files=self.full_path_concat(full_path_array) 

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
                Full_files_log=self.pd_resize_n(Simulation,simulation_id,search_id='simulation_id')
                reduced_files=Full_files_log[0]
                No_found_simu_id=Full_files_log[1]
                # path gene
                full_path_array=self.full_path(root_path,reduced_files)
                Selected_property_files=self.full_path_concat(full_path_array) 

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

    
            
        Prop=Selected_property_files
        self.neq_pandas=Prop[1]
        self.scatter_pandas=Prop[2]
        self.tmatrix_pandas=Prop[3]
        self.Selected_property_files=Selected_property_files

        return Selected_property_files,return_file_log,No_found_geo_id,No_found_mesh_id,No_found_material_id,No_found_simu_id

    def full_properties_look(self,neq_pandas=[],scatter_pandas=[],tmatrix_pandas=[],input=False):
        # get the file path 
        if input==False:
            neq_pandas=self.neq_pandas 
            scatter_pandas=self.scatter_pandas 
            tmatrix_pandas=self.tmatrix_pandas 
        else:
            pass
        
        neq_file_path=neq_pandas['full_path'].reset_index(drop=True)
        neq_len=len(neq_file_path)
        
        neq_simu_id=neq_pandas['simulation_id'].reset_index(drop=True)

        
        scatter_file_path=scatter_pandas['full_path'].reset_index(drop=True)
        scatter_len=len(scatter_file_path)

        scatter_simu_id=scatter_pandas['simulation_id'].reset_index(drop=True)

        
        tmatrix_file_path=tmatrix_pandas['full_path'].reset_index(drop=True)
        tmatrix_len=len(tmatrix_file_path)

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

            if  tmatrix_len==1:

                Tmatrix_properties=pd.DataFrame(Tmatrix_properties)
            else:
                Tmatrix_properties=pd.concat(Tmatrix_properties)
 
        else:
            Tmatrix_properties=pd.DataFrame()

        self.NEQ_properties=NEQ_properties
        self.Scatter_properties=Scatter_properties
        self.Tmatrix_properties=Tmatrix_properties
        return NEQ_properties, Scatter_properties, Tmatrix_properties

    def properties_look(self,properties,
                    NEQ_properties=[],
                    Scatter_properties=[],
                    Tmatrix_properties=[],group=False,input=False):
        '''
        Given a list of pandas arrays that contains full properties (the function: full_properties_look),
        we select the target properties from there. 
        (input: properties (neq_list,scatter_list,tmatirx_list), 2D list
                group: group by simulation if True        )
        
        
        MAHMOUD:  
            Please add description to this method here

            I suggest that this method should take only a single list of desired properties and return a single pandas DataFrame.
            If the properties is not found in a resultfile, then it is fine to ignore this resultfile

            I suggest also to have here an option or a default variable that can let the user decide whether to cluster all frequency- and/or angle-dependent properties in a sub DataFrame.
            For example, the final returned DataFrame should show the simulation_id only a single time, and at the row of the simulation_id, we need a sub DataFrame that has frequency- or angle-dependent properties,
            like [freq,scattering_cross_section]
        '''
        neq_properties=properties[0]
        scatter_properties=properties[1]
        tmatrix_properties=properties[2]
        
        # check is empty 
        if input==True:
            pass
        else:
            NEQ_properties=self.NEQ_properties
            Scatter_properties=self.Scatter_properties
            Tmatrix_properties=self.Tmatrix_properties
        neq_len=len(NEQ_properties)
        scatter_len=len(Scatter_properties)
        tmatrix_len=len(Tmatrix_properties)

        try:
            if neq_len!=0:
                # append simu id 
                neq_properties+=['simulation_id']
                Target_NEQ=NEQ_properties[neq_properties]
                #Target_NEQ.set_index(['simulation_id'],inplace=True)
            else:
                Target_NEQ=pd.DataFrame()
        except:
            Target_NEQ=['wrong properties name']

        try:
            if scatter_len!=0:
                # append simu id 
                
                scatter_properties+=['simulation_id']
        
                Target_Scatter=Scatter_properties[scatter_properties]
                #Target_Scatter.set_index(['simulation_id'],inplace=True)
            else:
                Target_Scatter=pd.DataFrame()
        except:
            Target_Scatter=['wrong properties name']

        try:
            if tmatrix_len!=0:
                # append simu id 
                tmatrix_properties+=['simulation_id']
                Target_Tmatrix=Tmatrix_properties[tmatrix_properties]
                #Target_Tmatrix.set_index(['simulation_id'],inplace=True)
            else:
                Target_Tmatrix=pd.DataFrame()
        except:
            Target_Tmatrix=['wrong properties name']
        if group:
            try:
                Target_NEQ=Target_NEQ.stack()
            except:
                pass
            try:
                Target_Scatter=Target_Scatter.stack()
            except:
                pass            
            try:
                Target_Tmatrix=Target_Tmatrix.stack()
            except:
                pass              

        return Target_NEQ, Target_Scatter, Target_Tmatrix
        

    def group_by(self,column,Check_list=['simulation_type', 'mesh_id', 'material_surr', 'material_particle', 'tolerance', 'RelOrient_id']):
       # try:
            # check list in 
        if column in Check_list:
            Simulation=self.Simulation
            Simu_path=self.full_path(self.root_path,Simulation)

            # make the set without the target search 
            Simu_path_drop=Simu_path.drop([column,'simulation_id','simulation_result_filepath','full_path','material_id'], axis=1).astype(str)
            # making new file that gives 
            All_rest_prop=Simu_path_drop.apply(lambda x: ''.join(x), axis=1)
            # create new concat frame
            New_concat_frame=pd.concat([All_rest_prop,Simu_path[column]],axis=1)
            # rename the column 
            New_concat_frame.columns = ['id',column ]
            # Adding simu id 
            New_concat_frame['simulation_id']=Simu_path['simulation_id']
            # check the counts for the unique group 
            New_unique_id_count=New_concat_frame.groupby('id')[column].count()
            unique_id=New_concat_frame['id'].unique()
            # get empty container to append each pandas array based on different column properties 
            column_container=[]
            for i in range(0,len(unique_id)):
                check_exist=New_concat_frame['id']==unique_id[i]
                simu_id_with_part_id=New_concat_frame['simulation_id'][check_exist]
                column_container.append(simu_id_with_part_id)

            self.unique_id=unique_id
            self.column_container=column_container

            return unique_id,column_container,New_unique_id_count
       # except:
      #      pass

    def get_prop_group_by(self,prop=[]):
        check_unique_contain=len(self.unique_id)
        Simu_container=self.column_container
        # open the file
        total_prop_list=[]
        #see_list=[]
        unique_id_others=self.unique_id
        
        for i in range(0,check_unique_contain):
            #current_see_list_for_one_unique_id=pd.Series(unique_id_others[i],index=['envi-id'])
            current_simu_id_list=Simu_container[i].reset_index(drop=True)
            # get len for simu_id
            current_simu_len=len(current_simu_id_list)
            current_simu_full_prop_list=[]
            for j in range(0,current_simu_len):
                current_simulation_id=current_simu_id_list[j]
                current_file_array=self.file_look(simulation_id=[current_simulation_id])
                # grab 
                neq_pandas=current_file_array[0][1]
                scatter_pandas=current_file_array[0][2]
                tmatrix_pandas=current_file_array[0][3]
                
                Full_prop=self.full_properties_look(neq_pandas=neq_pandas,scatter_pandas=scatter_pandas,tmatrix_pandas=tmatrix_pandas,input=True)
                
                if prop:
                    
                    Full_prop=self.properties_look(prop,Full_prop[0],Full_prop[1],Full_prop[2],input=True)
                    
                else:
                    pass
        
                current_nest_id_list=[Full_prop[0],Full_prop[1],Full_prop[2]]
                current_simu_full_prop_list.append(current_nest_id_list)
            total_prop_list.append(current_simu_full_prop_list)
        self.total_prop_list=total_prop_list
        return total_prop_list

            

    def data_see(self,envi_id,simu_type,prop):

        envi_id_loc=self.unique_id==envi_id
        envi_id_loc=[i for i, x in enumerate(envi_id_loc) if x]
        working_envi=self.total_prop_list[envi_id_loc[0]]
        print(len(working_envi))
        # type checking 
        if simu_type=='neq':
            type_number=0
        elif simu_type=='scatter':
            type_number=1
        elif simu_type=='tmatrix':
            type_number=2
        else:
            pass
        # plot
        fig = plt.figure()
        fig.set_figheight(15)
        fig.set_figwidth(15)

        ax = fig.add_subplot(111)
       # print(working_envi)
        for i in range(0,len(working_envi)):
            current_working_envi=working_envi[i][type_number].reset_index()
            print(current_working_envi)
            
            #x=current_working_envi[prop[0] ]
            #y=current_working_envi[prop[1] ]
            
            #ax.scatter(y,x,label=i)
            #ax.set_xlabel=prop[0]
            #ax.set_ylabel=prop[1]





        




# %%
