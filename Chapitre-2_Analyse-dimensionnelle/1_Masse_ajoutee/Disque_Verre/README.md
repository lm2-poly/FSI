# PyANSYS

## Folder structure

```images``` &rarr; Images necessary for ```added_mass.ipynb```

```saving_folder/disc_param_study``` &rarr; Folder for the results to be read in ```added_mass.ipynb```

```template-apdl-disc``` &rarr; Template input ```.mac``` files to launch simulations on Ansys for the discs

```template-apdl-glass``` &rarr; Template input ```.mac``` files to launch simulations on Ansys for the glass

```added_mass.ipynb``` &rarr; Jupyter Notebook for the post-treatment of results for the discs and glass to compute their added mass

```analyse_result_ansys``` &rarr; Storage and analysis of the results in adapted classes

```extract_result_ansys``` &rarr; Extraction of the results from a solved Ansys model

```launch_ansys_simulation``` &rarr; Code to launch a parametric simulation from a template input ```.mac``` file

```main_disc_analysis``` &rarr; Example of a parametric study for the disc

```main_glass_analysis``` &rarr; Example of a parametric study for the glass

## Requirements

To execute the code, it is needed to have Ansys APDL installed locally. The necessary libraries pyAnsys to communicate between Ansys and Python are given in ```requirements.txt``` in the home page, and the procedure to install them is available at the end of ```readme.md```.

## Objectives

This code allows the user to communicate with Ansys by the use of his own python code with the pyAnsys modules. 

One of the main goals of this code is launching, storing and analysing parametrical studies.

## How to generate a parametric study

Suppose that we have a ```.mac``` file working on APDL and we would like to modify it at each step to obtain a parametric study varying some parameters.

First we write a template of the ```.mac``` file that we want to modify. Then, we write a model class describing the model, one of the methods is writing the ```.mac``` file using the template. 

Then we set up the parametric study by creating a list of model objects containing the variation of the parameters we want. After this, we create the object "Parametric_study" with the list of models, it will generate automatically a list of results computed by Ansys under the "Analyse_result" object form. 

Finally you can analyse all the results you want from this "Parametric_study_object".
