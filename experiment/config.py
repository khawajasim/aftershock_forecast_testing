#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 15:15:28 2023

@author: khawajasim
"""


#1 - Run experiment with earthquakes in negative stress
introduce_FN = True 

#2 - Aggregated stress is already calculated from file "Data/calculated_stress.pkl". 
#    We can set it True to run it again
aggregation = False 

#3 - Save the results
save_results = True

#4 - Set to "True", if want to run simulation based on \DltaCFS along master fault.
#    Set to "False", if want to run simulation based on \DeltaCFS along optimally oriented fault.
stress_MAS= False   #True for MAS and False for OOP      

#5 - Model type to generate simulated catalogs.
#    Set 'model_type' = 'coulomb', to use Coulomb model for catalog generation.
#    Set 'model_type' = 'ref', to use Reference model for catalog generation. 
model_for_sim =  'coulomb' #'ref' 

#6 - Number of earthquakes in the simulated catalogs to run the analysis.
#Note: It can be any number for analysis, whether 200,300,.....2000.
total_eqs = 700 #2485  # #

#7 - Number of simulations to generate 
num_sim = 100


#8 - Grids to run the experiment for Chichi and Landers
grids = ['_single_res_grid_zoom=12','_single_res_grid_zoom=13','_single_res_grid_zoom=14','_multi_res_grid'] #
