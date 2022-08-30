#gestion des fichiers :
import pathlib
import os
#création de commandes terminal :
import argparse

#bibliothèque numpy :
import numpy as np

#importation de la bibliothèque créée :
import extract_result_ansys as ERA
import analyse_result_ansys as ARA
import launch_ansys_simulation as LAS

#Pour copier des variables :
from copy import copy
from copy import deepcopy

#pour faire de jolis print :
import pprint as pp

import matplotlib.pyplot as plt


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

args = parser.parse_args()

print("arguments of the parser :")
pp.pprint(vars(args))

n_points = args.point_number

#------------------------- The file paths --------------------------------
#-------------------------------------------------------------------------
print("------------init file paths--------------")
print("-----------------------------------------")
working_folder = pathlib.Path(os.getcwd())
param_folder = working_folder / "parametric-study"
param_folder.mkdir(parents = True, exist_ok = True)
saving_folder = working_folder / "saving_folder/glass_param_study"
saving_folder_no_water = working_folder / "saving_folder/glass_param_study_no_water"


model_name = {NOWATER:"cylinder_no_water.mac", WATER:"cylinder.mac"}

db_name = {NOWATER:"cylinder_no_water", WATER:"cylinder"}

template_file_path = {NOWATER:working_folder / "template-apdl-glass/cylinder_no_water.mac",
                      WATER:working_folder / "template-apdl-glass/cylinder.mac"}

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
average_radius = "average_radius"
width_structure = "width_structure"
outer_water_radius = "outer_water_radius"
radial_division = "radial_division"
angular_division = "angular_division"
DB_NAME = "db_name"
tol="tol"

if __name__ == "__main__":
    if args.main_case == "parametric_study":
        print("----beginning of parametric study -------")
        print("-----------------------------------------")
#------------------------- Variables init --------------------------------
#-------------------------------------------------------------------------
        print("---- initialisation of the variables-----")
        print("-----------------------------------------")
        #variables with water :
        variables = {}
        variables[average_radius] = 1.0
        variables[width_structure] = 0.01
        variables[outer_water_radius] = 2.0
        variables[radial_division] = 10 
        variables[angular_division] = 100
        variables[DB_NAME] = "'" + db_name[WATER] + "'"

        #variables without water :
        variables_no_water = variables.copy()
        variables_no_water[DB_NAME] = "'" + db_name[NOWATER] + "'"
        variables = {NOWATER:variables_no_water, WATER:variables}

        #other variables :
        fluid_element = "FLUID29"
        solid_element = "PLANE182"
        sound_celerity = 8000
        structural_density = 2500
        young_modulus = 70e9
        poisson = 0.22
        fluid_densitys = np.linspace(0,2000,n_points)

#------------------------- Result propreties -----------------------------
#-------------------------------------------------------------------------
        print("---- initialisation of the results propreties----")
        print("-------------------------------------------------")
        #result propreties with water :
        n_modes = 12
        modes_to_extract = [i for i in range(n_modes)]
        n_ddl = 3 #X, Y, pressure
        disp_ddls = [0,1]
        press_ddls = [2]
        fluid_comp = FLUID_COMP
        solid_comp = SOLID_COMP
        result_propreties = ARA.Result_propreties(modes_to_extract, n_ddl, solid_comp, fluid_comp, 
                                                disp_ddls, press_ddls)

        #result propreties without water :
        result_propreties_0 = deepcopy(result_propreties)
        n_modes = 10
        modes_to_extract = [i for i in range(n_modes)]
        result_propreties_0.n_ddl = 2
        result_propreties_0.press_ddls = None
        result_propreties_0.modes_to_extract = modes_to_extract
        result_propreties_0.extract_matrix = True
        result_propreties = {NOWATER:result_propreties_0, WATER:result_propreties}

#---------------------- Setting of models list ---------------------------
#-------------------------------------------------------------------------            
        print("---------- Setting of models list ---------")
        print("-------------------------------------------")
        models = []
        water = NOWATER
        for i, fluid_density in enumerate(fluid_densitys):
            if i > 0:
                water = WATER
            model = LAS.Glass_model(model_name[water], db_name[water], 
                        template_file_path[water], variables[water], fluid_element,
                        solid_element, fluid_density, sound_celerity,
                        structural_density, young_modulus, poisson, 
                        result_propreties[water])
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

if args.main_case == "no_water_extraction":
#------------------------- Variables init --------------------------------
#-------------------------------------------------------------------------
        print("---- initialisation of the variables-----")
        print("-----------------------------------------")
        #Variables into dictionnary :
        variables = {}
        variables[average_radius] = 1.0
        variables[width_structure] = 0.1
        variables[outer_water_radius] = 2.0
        variables[radial_division] = 10 
        variables[angular_division] = 10
        variables_no_water = variables.copy()
        variables_no_water[DB_NAME] = "'" + db_name[NOWATER] + "'"

        #other variables :
        fluid_element = "FLUID29"
        solid_element = "PLANE182"
        sound_celerity = 1430
        structural_density = 2500
        young_modulus = 70e9
        poisson = 0.22

#------------------------- Result propreties -----------------------------
#-------------------------------------------------------------------------
        print("---- initialisation of the results propreties----")
        print("-------------------------------------------------")
        #result propreties with water :
        n_modes = 20
        modes_to_extract = [i for i in range(n_modes)]
        n_ddl = 2 #X, Y
        disp_ddls = [0,1]
        press_ddls = None
        solid_comp = SOLID_COMP
        fluid_comp = None
        fluid_density = None
        result_propreties = ARA.Result_propreties(modes_to_extract, n_ddl, solid_comp, fluid_comp, 
                                                disp_ddls, press_ddls, extract_matrix=True)

        water = NOWATER
        print(model_name[water], db_name[water], result_file)
        model = LAS.Glass_model(model_name[water], db_name[water], 
                        template_file_path[water], variables, fluid_element,
                        solid_element, fluid_density, sound_celerity,
                        structural_density, young_modulus, poisson, 
                        result_propreties)
        
        models = [model]
#---------------------- The parametric study ---------------------------
#----------------------------------------------------------------------- 
        print("------------- Parametric study ------------")
        print("-------------------------------------------")
        parametric_study = LAS.Parametric_study(models)

#---------------------- Saving of the model ---------------------------
#----------------------------------------------------------------------
        print("------------ Saving the study :------------")
        print("-------------------------------------------")
        ERA.save_object(saving_folder_no_water, parametric_study)


if args.main_case == "load_no_water_study":
#------------------------- Variables init --------------------------------
#-------------------------------------------------------------------------
        print("--------- Loading of the study ---------")
        print("-----------------------------------------")
        parametric_study = ERA.load_object(saving_folder_no_water)
        parametric_study:LAS.Parametric_study

        results = parametric_study.results[0]
        results:ARA.Modal_result_w_matrixs
        K = results.K
        M = results.M
        # print(K)
        # print(M)
        # fig, (ax1, ax2) = plt.subplots(1, 2)
        # fig.suptitle("K and M Matrix profiles")
        # ax1.spy(K, markersize=0.01)
        # ax1.set_title("K Matrix")
        # ax2.spy(M, markersize=0.01)
        # ax2.set_title("M Matrix")
        # plt.show(block=True)
        i = 2
        eigen_vec = results.get_eigen_vec(i)
        modal_stiff = results.modal_stiff(i)
        modal_mass = results.modal_mass(i)
        ev = results.eigen_val(i)
        print(eigen_vec)
        print(modal_mass)
        print(modal_stiff)
        print(modal_stiff - ev * modal_mass)

if args.main_case == "load_full_study":
#------------------------- Variables init --------------------------------
#-------------------------------------------------------------------------
        print("--------- Loading of the study ---------")
        print("-----------------------------------------")
        parametric_study = ERA.load_object(saving_folder)
        parametric_study:LAS.Parametric_study

        print("--------- Extracting the results ---------")
        print("-----------------------------------------")
        results = parametric_study.results
        structural_result = results[0]
        structural_result:ARA.Modal_result_w_matrixs
        structural_freqs, idx_str_freq = structural_result.get_frequencies(-1, True)
        nb_modes = 5
        structural_freqs, idx_str_freq = structural_freqs[:nb_modes], idx_str_freq[:nb_modes]
        print(structural_freqs)
        results_fluid = results[1:]
        last_freqs = copy(structural_freqs)
        print("--------- Computing of the added mass ---------")
        print("-----------------------------------------")
        added_mass = np.zeros((len(results_fluid) + 1, len(structural_freqs)))
        frequency_selected = np.zeros((len(results_fluid) + 1, len(structural_freqs)))
        frequency_selected[0,:] = structural_freqs
        scores = np.zeros((len(results_fluid), ))
        for idx_param, result_fluid in enumerate(results_fluid):
                result_fluid:ARA.Modal_result
                approxs_freq = np.zeros((len(structural_freqs),))
                for idx_freq, struc_freq in enumerate(structural_freqs):
                        idx_struc_freq = idx_str_freq[idx_freq]
                        if idx_param < 2:
                                approx_freq = last_freqs[idx_freq]
  
                        else:
                                freq_1 = frequency_selected[idx_param - 1, idx_freq]
                                freq_0 = frequency_selected[idx_param - 2, idx_freq]
                                approx_freq = 2 * freq_1 - freq_0
                        approxs_freq[idx_freq] = approx_freq
                                      
                freqs, struct_ids, score = result_fluid.associate_frequencies(approx_freqs=approxs_freq, 
                score = True, family_id=True)
                scores[idx_param] = score
                last_freqs = freqs    
                print(struct_ids)
                for idx, freq in enumerate(freqs):
                        idx_struc_freq = idx_str_freq[struct_ids[idx]]
                        
                        added_mass[idx_param + 1, idx] = \
                                structural_result.compute_added_mass(freq, 
                                idx_struc_freq)
                frequency_selected[idx_param + 1, :] = freqs
        structural_result.supp_matrixs()
        # print(frequency_selected)
        print(added_mass)
        print(frequency_selected)
        print(scores)
        densitys = parametric_study.variable_list(LAS.DENS)
        plt.plot(densitys, added_mass)
        plt.show()
        print("------------ Saving the study :------------")
        print("-------------------------------------------")
        ERA.save_object(saving_folder, parametric_study)

                        



        # idx_mode_to_analyse = [0, 3, 5, 6]
        # idx_mode_to_analyse_with_water = [0, 3, 7, 9]
        # results = results[1:]
        # for i, idx_mode in enumerate(idx_mode_to_analyse):
        #         added_mass = []
        #         for result in results:
        #                 result:ARA.Modal_result
        #                 frequency = result.frequencies[idx_mode_to_analyse_with_water[i]]
        #                 added_mass_ = result_without_water.compute_added_mass(
        #                         frequency, idx_mode)
        #                 added_mass.append(added_mass_)

        #         print("---------------- Plotting ---------------")
        #         print("-----------------------------------------")
        #         densitys = parametric_study.variable_list(LAS.DENS)
        #         densitys = densitys[1:]
        #         plt.plot(densitys, added_mass)
        # plt.show()