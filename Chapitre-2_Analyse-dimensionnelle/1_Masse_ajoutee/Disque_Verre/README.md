# PyANSYS

## Requirements

To use the code, it is needed to have Ansys APDL installed locally. The necessary libraries pyAnsys to communicate between Ansys and Python are given in ```requirements.txt``` in the home page, and the procedure to install them is available at the end of ```readme.md```.

## Objectives

This code allows the user to communicate with Ansys by the use of his own python code with the pyAnsys modules. 

One of the main goals of this code is launching, storing and analysing parametrical studies.

## Description of the code structure

The code is divided in 3 main modules:
* extract_result_ansys: to extract the results from a solved ansys model
* analyse_result_ansys: to store and analyse the results in adapted classes
* launch_ansys_simulation: to launch a parametric simulation from a template input ```.mac```

## How to generate a parametric study

Suppose that we have a ```.mac``` file working on APDL and we would like to modify it at each step to obtain a parametric study varying some parameters.

First we write a template of the ```.mac``` file that we want to modify. Then, we write a model class describing the model, one of the methods is writing the ```.mac``` file using the template. 

Then we set up the parametric study by creating a list of model objects containing the variation of the parameters we want. After this, we create the object "Parametric_study" with the list of models, it will generate automatically a list of results computed by Ansys under the "Analyse_result" object form. 

Finally you can analyse all the results you want from this "Parametric_study_object".

### Examples of a parametric study available

We gives two examples of a parametric study processed using this method, all the steps are described in the 2 files:
* main_disc_analysis
* main_glass_analysis 

Also the results of these parametric studies are analysed with the aim of computing the added mass of different systems.
