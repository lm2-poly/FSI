#gestion des fichiers :
import pathlib
import os

#création de commandes terminal :
import argparse

#bibliothèque numpy et matplotlib :
import numpy as np
import matplotlib.pyplot as plt

#importation de la bibliothèque créée :
import extract_result_ansys as ERA
import analyse_result_ansys as ARA
import launch_ansys_simulation as LAS

#Pour copier des variables :
from copy import deepcopy

#pour faire de jolis print :
import pprint as pp


#------------------------- Global constants ------------------------------
#-------------------------------------------------------------------------
print("----------init global constants----------")
print("-----------------------------------------")
FREQ = "freq"
COMP = "comp"
NODAL = "nodal-result"
POS = "nodal-position"
ID_INDEX = "id-to-index"
FREQ_INDEX = "freq-index"
SOLID_COMP = "SOLID_NODES"
FLUID_COMP = "FLUID_NODES"
WATER = 1
NOWATER = 0

#------------------------- Argument parser -------------------------------
#-------------------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument('--main_case')
parser.add_argument('--point_number', type = int, default=20)
parser.add_argument('--is_unconfined', action='store_true', default=False)

args = parser.parse_args()

print("arguments of the parser :")
pp.pprint(vars(args))

if args.is_unconfined:
    analysis_type = "disc_unconfined"
else:
    analysis_type = "disc_confined"

n_points = args.point_number

#------------------------- The file paths --------------------------------
#-------------------------------------------------------------------------
print("------------init file paths--------------")
print("-----------------------------------------")
working_folder = pathlib.Path(os.getcwd())
param_folder = working_folder / "parametric-study"
param_folder.mkdir(parents = True, exist_ok = True)
print(analysis_type)
saving_folder = working_folder / "saving_folder/disc_param_study" / analysis_type

model_name = {NOWATER:"disc_no_water.mac", WATER:"disc.mac"}

db_name = {NOWATER:"disc_no_water", WATER:"disc"}

template_file_path = {NOWATER:working_folder / ("template-apdl-disc" +  
                "/" + analysis_type + "/" + "/disc_no_water.mac"),
                      WATER:working_folder / ("template-apdl-disc" +  
                "/" + analysis_type + "/" + "/disc.mac")}

result_file = {NOWATER:db_name[NOWATER] + ".rst",
               WATER: db_name[WATER] + ".rst"}
result_file = {NOWATER: param_folder / result_file[NOWATER],
               WATER: param_folder / result_file[WATER]}         

data_file = {NOWATER:db_name[NOWATER] + ".dat",
               WATER: db_name[WATER] + ".dat"}
data_file = {NOWATER: param_folder / data_file[NOWATER],
               WATER: param_folder / data_file[WATER]}        

#------------------------- Variables name --------------------------------
#-------------------------------------------------------------------------
print("----------init variables name------------")
print("-----------------------------------------")

inner_radius = "inner_radius"
outer_radius = "outer_radius"
spring_length="spring_length"
width_structure = "width_structure"
outer_water_radius = "outer_water_radius"
radial_division = "radial_division"
angular_division = "angular_division"
DB_NAME = "db_name"
tol="tol"
r1 = "real_constant_1" #spring stiffness
r2 = "real_constant_2" #mass value
r3 = "real_constant_3" #reference pressure for acoustic_element
r4 = "real_constant_4" #RAD, X0, Y0 for FLUID129


if __name__ == "__main__":
    if args.main_case == "parametric_study":
        print("----" + "case : " + args.main_case + "----")
        print("----beginning of parametric study -------")
        print("-----------------------------------------")
#------------------------- Variables init --------------------------------
#-------------------------------------------------------------------------
        print("---- initialisation of the variables-----")
        print("-----------------------------------------")
        #variables with water :
        variables = {}
        variables[inner_radius] = .2
        variables[outer_radius] = 0.4
        variables[width_structure] = 0.01
        variables[spring_length] = 0.6
        variables[radial_division] = 10 
        variables[angular_division] = 10
        variables[tol]=0.0015
        variables[DB_NAME] = "'" + db_name[WATER] + "'"

        #variables without water :
        variables_no_water = variables.copy()
        variables_no_water[DB_NAME] = "'" + db_name[NOWATER] + "'"
        variables = {NOWATER:variables_no_water, WATER:variables}

        #other variables :
        mass_element = "MASS21"
        spring_element = "COMBIN14"
        fluid_element = "FLUID29"
        fluid_infinite_element = "FLUID129"
        solid_element = "PLANE182"
        sound_celerity = 1600
        structural_sound_celerity = 5000
        structural_density = 2700
        young_modulus = 69000000000
        poisson = 0.346
        fluid_densitys = np.linspace(0,2000,n_points)
        real_constants={}
        real_constants[r1]=40000 #spring_stiffness
        real_constants[r2]= 200 #mass
        real_constants[r3]=2e-5 #reference pressure for acoustic_element
        real_constants[r4]=.4

#------------------------- Result propreties -----------------------------
#-------------------------------------------------------------------------
        print("---- initialisation of the results propreties----")
        print("-------------------------------------------------")
        #result propreties with water :
        n_modes = 4
        modes_to_extract = [i for i in range(n_modes)]
        if args.is_unconfined:
            n_ddl = 9
        else:
            n_ddl = 7
        disp_ddls = [0,1,2]
        press_ddls = [3]
        fluid_comp = FLUID_COMP
        solid_comp = SOLID_COMP
        result_propreties = ARA.Result_propreties(modes_to_extract, n_ddl, solid_comp, fluid_comp, 
                                                disp_ddls, press_ddls)

        #result propreties without water :
        result_propreties_0 = deepcopy(result_propreties)
        n_modes = 1
        modes_to_extract = [i for i in range(n_modes)]
        result_propreties_0.n_ddl = 7
        result_propreties_0.press_ddls = None
        result_propreties_0.modes_to_extract = modes_to_extract
        result_propreties = {NOWATER:result_propreties_0, WATER:result_propreties}

#---------------------- Setting of models list ---------------------------
#-------------------------------------------------------------------------            
        print("---------- Setting of models list ---------")
        print("-------------------------------------------")
        models = []
        water = NOWATER
        for i, fluid_density  in enumerate(fluid_densitys):
            if i > 0:
                water = WATER
            model = LAS.Disc_model(model_name[water], db_name[water], 
        template_file_path[water], variables[water], fluid_element,
        solid_element, fluid_density, sound_celerity, structural_density,
        young_modulus, poisson, result_propreties[water], real_constants, 
        mass_element, spring_element, fluid_infinite_element, 
        structural_sound_celerity)
            models.append(model)

#---------------------- The parametric study ---------------------------
#----------------------------------------------------------------------- 
        print("------------- Parametric study ------------")
        print("-------------------------------------------")
        parametric_study = LAS.Parametric_study(models)

#---------------------- Saving of the model ---------------------------
#----------------------------------------------------------------------
        print("------------ Saving the study :------------")
        print("-------------------------------------------")
        ERA.save_object(saving_folder, parametric_study)
        print("saved in :", saving_folder)
        
    if args.main_case == "analyse-result-disc":
        print("----" + "case : " + args.main_case + " ----")
        parametric_result:LAS.Parametric_study = ERA.load_object(saving_folder)
        n_models = parametric_result.n_models()
        models = parametric_result.models
        results = parametric_result.results
        #computation of the added mass:
        added_mass = np.zeros((n_models,))
        for i, (model, result) in enumerate(zip(models, results)):
            print("case",i)
            #declare the model and the result to help 
            #to write the code (autocompletion)
            model:LAS.Disc_model #declaration des variables
            result:ARA.Modal_result
            frequencies, ids = result.get_frequencies(decimal = 4, ids_ = True)
            if i == 0:
                dry_mode_frequency = frequencies[0]
                frequency=frequencies[0]
            else:
                frequency = frequencies[1]
            print("freq :", frequency)
            spring_stiffness = model.real_constants[r1]
            added_mass[i] = (spring_stiffness / (4 * np.pi ** 2) 
            * (1 / frequency ** 2 - 1 / dry_mode_frequency ** 2))
        print("added mass :",added_mass)
        fluid_densitys = parametric_result.variable_list("density")
        plt.figure()
        title="frequency :"+str(dry_mode_frequency)+' Hz'
        plt.plot(fluid_densitys,added_mass)
        plt.title(title)
        plt.show()
        
        models = parametric_result.models



        
        
        