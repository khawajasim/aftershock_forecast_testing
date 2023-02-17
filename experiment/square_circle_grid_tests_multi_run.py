#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 16:08:12 2022

@author: khawajasim
"""
import scipy
import numpy
import numpy as np
import pandas
import matplotlib.pyplot as plt
from shapely.geometry import box, shape, mapping, Point, Polygon
import json
from grid_operations import bounds_to_polygons, create_square_grid_bounds, create_circular_grid, forecast_aggregation
from evaluations import calc_ROC_MCCF1, calc_mccf1_metric
from utils import write_geojson_feature, generate_gridded_eqs
from labellines import labelLines

introduce_FN = True
aggregation = False
save_results = True #False
stress_MAS = True  #True for MAS and False for OOP
model_type =  'coulomb' #'ref' 
num_sim = 100
total_eqs = 400 #2485  # 



if aggregation:
    #---------------Start of code
    events = pandas.read_pickle('coulomb_forecast/Event_4.pkl')
    stress_data = events.to_numpy()
    depth = 7 #5 if Event_3
    stress_data = stress_data[stress_data[:,2] == depth]
    
    #Rescaling the degrees to the km distance --- -150 to 150
    stress_data[:,0] = numpy.round(stress_data[:,0]/max(stress_data[:,0]) * 150)
    stress_data[:,1] = numpy.round(stress_data[:,1]/max(stress_data[:,1]) * 150)
    
    #Filter within 100km x 100km
    within_range = numpy.logical_and(numpy.logical_and(stress_data[:,0]>=-100, stress_data[:,1]>=-100),
                      numpy.logical_and(stress_data[:,0]<100, stress_data[:,1]<100))
    
    stress_data= stress_data[within_range]
    
    #####-========------ Convert the forecast points into grid cells: MAKE A GEOJSON for QGIS----
    #Considering the given points as origin coordinates
    dh = numpy.diff(numpy.unique(stress_data[:,1]))[0]
    cell_bounds = numpy.column_stack((stress_data[:,:2], stress_data[:,:2]+dh))
    model_grid = bounds_to_polygons(cell_bounds)
    filename = 'square_circle/model_grid_mas_rescaledCoords'
#    write_geojson_feature(model_grid,stress_data[:,3],filename )
    ###--========-- Convert the forecast points into grid cells here .....
    if stress_MAS:
        stress_data = numpy.column_stack((stress_data[:,:3], stress_data[:,3])) 
    else:
        stress_data = numpy.column_stack((stress_data[:,:3], stress_data[:,4]))
    ####--------NEXT STEPS FROM HERE------ (10-10-2022)
    # 1 -  Generate distance based stress
    
    #Choose location data and OOP stress
#    stress_data = numpy.column_stack((stress_data[:,:3], stress_data[:,3])) #4
    
    # # Deleted the centeral point of the whole grid, to make it symetric for plotting with "plot_stress" function.
    # # And secondly to calculate the static stress from origin point (0,0). Because distance from origin point is zero, leading to inf stress.
    # stress_data = numpy.delete(stress_data, numpy.where(numpy.logical_or(stress_data[:,0]==0, stress_data[:,1]==0) == True), 0 ) 
    
    #data = numpy.loadtxt('stress_oop.csv', skiprows=1, delimiter=',')
    
     #-- Generate Static Stress
     # center = [0,0]
    coords = stress_data[:,:2]
    dist = []
    for coord in coords:
        dist.append(numpy.sqrt(coord[0]**2 + coord[1]**2))
        
    dist = numpy.array(dist) 
    c = 2 #Static stress constant
    static_stress = c* dist ** -2  #-3 (2 by Hardebeck)
    stress_data = numpy.column_stack((stress_data, static_stress))
    
    if numpy.isinf(max(static_stress)):
        static_stress[numpy.where(numpy.isinf(static_stress) == True)] = numpy.sort(static_stress)[-2]
#        stress_data = numpy.delete(stress_data, numpy.where(numpy.isinf(static_stress) == True),0)
#        cell_bounds = numpy.delete(cell_bounds, numpy.where(numpy.isinf(static_stress) == True),0)
        
    
    model_grid = bounds_to_polygons(cell_bounds)
    filename = 'forecast_data/model_grid_ref'
    write_geojson_feature(model_grid,stress_data[:,4],filename)


#-----SQUARE TEST GRID
square_bounds = create_square_grid_bounds([-100,-100], [100,100],40)
square_grid = bounds_to_polygons(square_bounds)
#---- CIRCLE TEST GRID
org = (0,0)
r_seg = 40
a_seg= 40
radius_max  = 100
circle_grid = create_circular_grid(radius_max, r_seg, a_seg, origin = org)


if stress_MAS:
        fname_square = 'forecast_data/square_grid_aggregated_stress_'+str(len(square_grid))+'_MAS.csv'
        fname_circle = 'forecast_data/circle_grid_aggregated_stress_'+str(len(circle_grid))+'_MAS.csv'
else:
        fname_square = 'forecast_data/square_grid_aggregated_stress_'+str(len(square_grid))+'_OOP.csv'
        fname_circle = 'forecast_data/circle_grid_aggregated_stress_'+str(len(circle_grid))+'_OOP.csv'


if aggregation:
    square_grid_oop = forecast_aggregation(square_grid, model_grid, stress_data[:,3])
    write_geojson_feature(square_grid,square_grid_oop,'forecast_data/square_grid_mas_'+str(numpy.size(square_grid)))
    
    square_grid_ref = forecast_aggregation(square_grid, model_grid, stress_data[:,4])
    write_geojson_feature(square_grid, square_grid_ref, 'forecast_data/square_grid_ref_'+str(numpy.size(square_grid)))
    
    stress_data_square = numpy.column_stack((square_grid_oop, square_grid_ref))
    
    
    circle_grid_oop = forecast_aggregation(circle_grid, model_grid, stress_data[:,3])
    write_geojson_feature(circle_grid, circle_grid_oop,'forecast_data/circle_grid_mas_'+str(numpy.size(circle_grid)))
    
    circle_grid_ref = forecast_aggregation(circle_grid, model_grid, stress_data[:,4])
    write_geojson_feature(circle_grid, circle_grid_ref,'forecast_data/circle_grid_ref_'+str(numpy.size(circle_grid)))
    
    stress_data_circle = numpy.column_stack((circle_grid_oop, circle_grid_ref))
    numpy.savetxt(fname_square, 
                  stress_data_square, delimiter=',', header='MAS,Ref', comments='')
    numpy.savetxt(fname_circle, 
                  stress_data_circle, delimiter=',', header='MAS,Ref', comments='')

#------ Locate earthquakes to grid ---
# numpy.savetxt('stress_data_square_test_grid_'+str(len(stress_data_square))+'.csv', stress_data_square, delimiter=',', header='oop,ref', comments='')
# numpy.savetxt('stress_data_circle_test_grid_'+str(len(stress_data_circle))+'.csv', stress_data_circle, delimiter=',', header='oop,ref', comments='')

stress_data_square = numpy.loadtxt(fname_square , delimiter=',', skiprows = 1)
stress_data_circle = numpy.loadtxt(fname_circle, delimiter=',', skiprows =1)


#-- Generate synthetic catalog
#Use probability map of static stress to generate catalog. 
#But generate events in only those cells which have positive OOP.
#The number obtained by searching http://www.isc.ac.uk/iscbulletin/search/catalogue/
#for 2 months after Amaratice earthquake

FNs = [0,1, 2,3,4,5,6,7]

roc_oop_FNs_square = 1*FNs
roc_ref_FNs_square = 1*FNs
mcc_f1_oop_FNs_square = 1*FNs
mcc_f1_ref_FNs_square = 1*FNs

roc_oop_FNs_circ = 1*FNs
roc_ref_FNs_circ = 1*FNs
mcc_f1_oop_FNs_circ = 1*FNs
mcc_f1_ref_FNs_circ = 1*FNs

for i in range(num_sim):
    #---- FNs = 6
    roc_oop_square = []
    roc_ref_square = []
    mcc_f1_oop_square = []
    mcc_f1_ref_square = []
    
    roc_oop_circ = []
    roc_ref_circ = []
    mcc_f1_oop_circ = []
    mcc_f1_ref_circ = []
# ---------
#for grid in grid_types:
    #---- FNs = 6
    for fn in FNs:
        print('False Negative :',fn)
        square_grid_arranged, stress_data_square_arranged = generate_gridded_eqs(square_grid, stress_data_square, 
                                                                                 model_type=model_type, Num_eqs=[total_eqs-fn, fn])
        #write_geojson_feature(square_grid_arranged.tolist(), stress_data_square_arranged[:,2],
        #'results_new/square_grid_eqs_arranged_'+str(fn))
        #---
        circle_grid_arranged, stress_data_circle_arranged = generate_gridded_eqs(circle_grid, stress_data_circle, 
                                                                                 model_type=model_type, Num_eqs=[total_eqs-fn, fn])
        #write_geojson_feature(circle_grid_arranged.tolist(), stress_data_circle_arranged[:,2],
        #                       'results_new/circle_grid_eqs_arranged_'+str(fn))
    
        #----============== PLOTS for Square Grids.
        oop_stress = stress_data_square_arranged[:,0]
        static_stress = stress_data_square_arranged[:,1]
        earthquakes = stress_data_square_arranged[:,2]
        print('No. of Active cells for SQUARE Grid:', len(earthquakes[earthquakes>0]))
        print('----Active FNs cells in SQUARE Grid:', len(oop_stress[numpy.logical_and(earthquakes>0, oop_stress<0)]))
        
        #--- ------Plot ROC OOP stress
        oop_ROCdata, oop_MCC_F1 = calc_ROC_MCCF1(oop_stress, earthquakes)
        fpr = oop_ROCdata[:,0]
        tpr = oop_ROCdata[:,1]
        auc = numpy.round(abs(numpy.trapz(tpr, fpr)), 3)
         
        roc_oop_square.append(auc)
        f1 = oop_MCC_F1[:,0]
        mcc = oop_MCC_F1[:,1]
        mcc_f1 = numpy.round(calc_mccf1_metric(f1,mcc, min_dist=True),3)
         
        mcc_f1_oop_square.append(mcc_f1)
        static_ROCdata, static_MCC_F1 = calc_ROC_MCCF1(static_stress, earthquakes)
        fpr = static_ROCdata[:,0]
        tpr = static_ROCdata[:,1]
        auc = numpy.round(abs(numpy.trapz(tpr, fpr)), 3)
        roc_ref_square.append(auc)
              
        f1 = static_MCC_F1[:,0]
        mcc = static_MCC_F1[:,1]
        mcc_f1 = calc_mccf1_metric(f1,mcc, min_dist=True)
        mcc_f1_ref_square.append(mcc_f1)
            
            
        #----============== PLOTS for CIRCULAR Grids.
        oop_stress = stress_data_circle_arranged[:,0]
        static_stress = stress_data_circle_arranged[:,1]
        earthquakes = stress_data_circle_arranged[:,2]
        print('No. of Active cells for CIRCULAR Grid:', len(earthquakes[earthquakes>0]))
        print('----Active FNs cells in CIRCULAR Grid:', len(oop_stress[numpy.logical_and(earthquakes>0, oop_stress<0)]))
    
            
        #--- ------Plot ROC OOP stress
        oop_ROCdata, oop_MCC_F1 = calc_ROC_MCCF1(oop_stress, earthquakes)
        fpr = oop_ROCdata[:,0]
        tpr = oop_ROCdata[:,1]
        auc = numpy.round(abs(numpy.trapz(tpr, fpr)), 3)
        roc_oop_circ.append(auc)
         
        f1 = oop_MCC_F1[:,0]
        mcc = oop_MCC_F1[:,1]
        mcc_f1 = numpy.round(calc_mccf1_metric(f1,mcc, min_dist=True),3)
        mcc_f1_oop_circ.append(mcc_f1)
        static_ROCdata, static_MCC_F1 = calc_ROC_MCCF1(static_stress, earthquakes)
        fpr = static_ROCdata[:,0]
        tpr = static_ROCdata[:,1]
        auc = numpy.round(abs(numpy.trapz(tpr, fpr)), 3)
        roc_ref_circ.append(auc)
                
        f1 = static_MCC_F1[:,0]
        mcc = static_MCC_F1[:,1]
        mcc_f1 = numpy.round(calc_mccf1_metric(f1,mcc, min_dist=True),3)
        mcc_f1_ref_circ.append(mcc_f1)
        
    roc_oop_FNs_square = numpy.row_stack((roc_oop_FNs_square, roc_oop_square))
    mcc_f1_oop_FNs_square = numpy.row_stack((mcc_f1_oop_FNs_square, mcc_f1_oop_square))
    roc_ref_FNs_square = numpy.row_stack((roc_ref_FNs_square, roc_ref_square))
    mcc_f1_ref_FNs_square = numpy.row_stack((mcc_f1_ref_FNs_square, mcc_f1_ref_square))
    
    
    roc_oop_FNs_circ = numpy.row_stack((roc_oop_FNs_circ, roc_oop_circ))
    mcc_f1_oop_FNs_circ = numpy.row_stack((mcc_f1_oop_FNs_circ, mcc_f1_oop_circ))
    roc_ref_FNs_circ = numpy.row_stack((roc_ref_FNs_circ, roc_ref_circ))
    mcc_f1_ref_FNs_circ = numpy.row_stack((mcc_f1_ref_FNs_circ, mcc_f1_ref_circ))
        
      
# #------ Square ROC        
# fig1, axs1 = plt.subplots()
# fig1.set_size_inches(20, 15)
# # axs1.hist([roc_oop_FNs_square[1:,0],roc_oop_FNs_square[1:,1], roc_oop_FNs_square[1:,2], roc_oop_FNs_square[1:,3],roc_oop_FNs_square[1:,4],
# #            roc_oop_FNs_square[1:,5], roc_oop_FNs_square[1:,6]], 
# #          label=['FN=0', 'FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'],
# #          color=['b', 'g', 'r','y','m','c','k']) #, roc_oop_FNs[1:,7]]--- --

# # axs1.hist([roc_ref_FNs_square[1:,0],roc_ref_FNs_square[1:,1], roc_ref_FNs_square[1:,2], roc_ref_FNs_square[1:,3],roc_ref_FNs_square[1:,4],
# #            roc_ref_FNs_square[1:,5], roc_ref_FNs_square[1:,6]], 
# #          label=['FN=0','FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'],
# #          color=['b', 'g', 'r','y','m','c','k'], 
# #          histtype='step', linewidth=2)    # ---,'FN=7'---,
# axs1.hist(roc_oop_FNs_square[1:,0], label='Coulomb Model', color='r') #, roc_oop_FNs[1:,7]] 
# axs1.hist(roc_ref_FNs_square[1:,0], label='Reference Model', color='g', linewidth=2)    # ---,'FN=7'---,histtype='step', 

# axs1.legend()

# axs1.hist([roc_oop_FNs_square[1:,1], roc_oop_FNs_square[1:,2], roc_oop_FNs_square[1:,3],roc_oop_FNs_square[1:,4],
#            roc_oop_FNs_square[1:,5], roc_oop_FNs_square[1:,6]], 
#          label=['FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'], alpha=0.3, color=[ 'r','r','r','r','r','r'], 
#          linewidth=2, linestyle='--')    # ---,'FN=7'---,histtype='step', 



# axs1.set_title('Square Grid - AUC value: Reference Model vs Coulomb Model')
# axs1.set_xlabel('AUC value')

# #------Circular ROC

# fig2, axs2 = plt.subplots()
# fig2.set_size_inches(20, 15)

# # axs2.hist([roc_oop_FNs_circ[1:,0],roc_oop_FNs_circ[1:,1], roc_oop_FNs_circ[1:,2], roc_oop_FNs_circ[1:,3],roc_oop_FNs_circ[1:,4],
# #            roc_oop_FNs_circ[1:,5], roc_oop_FNs_circ[1:,6]], 
# #          label=['FN=0', 'FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'],
# #          color=['b', 'g', 'r','y','m','c','k']) #, roc_oop_FNs[1:,7]]--- -- 
# # axs2.legend()
# # axs2.hist([roc_ref_FNs_circ[1:,0],roc_ref_FNs_circ[1:,1], roc_ref_FNs_circ[1:,2], roc_ref_FNs_circ[1:,3],roc_ref_FNs_circ[1:,4],
# #            roc_ref_FNs_circ[1:,5], roc_ref_FNs_circ[1:,6]], 
# #          label=['FN=0','FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'],
# #          color=['b', 'g', 'r','y','m','c','k'], 
# #          histtype='step', linewidth=2)    # ---,'FN=7'---,

# axs2.hist(roc_oop_FNs_circ[1:,0], label='Coulomb Mobel', color='r') #, roc_oop_FNs[1:,7]]--- -- 
# axs2.hist(roc_ref_FNs_circ[1:,0], label='Reference Model', color='g')    # ---,'FN=7'---,histtype='step', , linewidth=2 

# axs2.legend()
# axs2.hist([roc_oop_FNs_circ[1:,1], roc_oop_FNs_circ[1:,2], roc_oop_FNs_circ[1:,3],roc_oop_FNs_circ[1:,4],
#            roc_oop_FNs_circ[1:,5], roc_oop_FNs_circ[1:,6]], 
#          label=['FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'], alpha=0.3, color=[ 'r','r','r','r','r','r'], 
#          linewidth=2, linestyle='--')    # ---,'FN=7'---,histtype='step', 

# axs2.set_title('Circle Grid - AUC value: Reference Model vs Coulomb Model')
# axs2.set_xlabel('AUC value')


# #----------SQUARE MCC-F1
# fig3, axs3 = plt.subplots()
# fig3.set_size_inches(20, 15)

# # axs3.hist([mcc_f1_oop_FNs_square[1:,0],mcc_f1_oop_FNs_square[1:,1], mcc_f1_oop_FNs_square[1:,2], mcc_f1_oop_FNs_square[1:,3],mcc_f1_oop_FNs_square[1:,4],
# #            mcc_f1_oop_FNs_square[1:,5], mcc_f1_oop_FNs_square[1:,6]], 
# #          label=['FN=0', 'FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'],
# #          color=['b', 'g', 'r','y','m','c','k']) #, roc_oop_FNs[1:,7]]--- -- 
# # axs3.legend()
# # axs3.hist([mcc_f1_ref_FNs_square[1:,0],mcc_f1_ref_FNs_square[1:,1], mcc_f1_ref_FNs_square[1:,2], mcc_f1_ref_FNs_square[1:,3],mcc_f1_ref_FNs_square[1:,4],
# #            mcc_f1_ref_FNs_square[1:,5], mcc_f1_ref_FNs_square[1:,6]], 
# #          label=['FN=0','FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'],
# #          color=['b', 'g', 'r','y','m','c','k'], 
# #          histtype='step', linewidth=2)    # ---,'FN=7'---,

# axs3.hist(mcc_f1_oop_FNs_square[1:,0], label='Coulomb Model', color='r') #, roc_oop_FNs[1:,7]]--- -- 
# axs3.hist(mcc_f1_ref_FNs_square[1:,0], label='Reference Model', color='g')    #histtype='step', , linewidth=2 
# axs3.legend()
# axs3.hist([mcc_f1_oop_FNs_square[1:,1], mcc_f1_oop_FNs_square[1:,2], mcc_f1_oop_FNs_square[1:,3],mcc_f1_oop_FNs_square[1:,4],
#             mcc_f1_oop_FNs_square[1:,5], mcc_f1_oop_FNs_square[1:,6]], 
#           label=['FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'], alpha=0.3, color=['r', 'r', 'r','r','r','r'],
#           linewidth=2, linestyle='--')    # ---,'FN=7'---,histtype='step', #, roc_oop_FNs[1:,7]]--- -- 


# axs3.set_title('Square Grid - MCC-F1 value: Reference Model vs Coulomb Model')
# axs3.set_xlabel('MCf1 value')
        

# #------
# #---------- CIRCLE MCC-F1
# fig4, axs4 = plt.subplots()
# fig4.set_size_inches(20, 15)

# # axs4.hist([mcc_f1_oop_FNs_circ[1:,0],mcc_f1_oop_FNs_circ[1:,1], mcc_f1_oop_FNs_circ[1:,2], mcc_f1_oop_FNs_circ[1:,3],mcc_f1_oop_FNs_circ[1:,4],
# #            mcc_f1_oop_FNs_circ[1:,5], mcc_f1_oop_FNs_circ[1:,6]], 
# #          label=['FN=0', 'FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'],
# #          color=['b', 'g', 'r','y','m','c','k']) #, roc_oop_FNs[1:,7]]--- -- 
# # axs4.legend()
# # axs4.hist([mcc_f1_ref_FNs_circ[1:,0],mcc_f1_ref_FNs_circ[1:,1], mcc_f1_ref_FNs_circ[1:,2], mcc_f1_ref_FNs_circ[1:,3],mcc_f1_ref_FNs_circ[1:,4],
# #            mcc_f1_ref_FNs_circ[1:,5], mcc_f1_ref_FNs_circ[1:,6]], 
# #          label=['FN=0','FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'],
# #          color=['b', 'g', 'r','y','m','c','k'], 
# #          histtype='step', linewidth=2)    # ---,'FN=7'---,

# axs4.hist(mcc_f1_oop_FNs_circ[1:,0], label='Coulomb Model', color='r') #, roc_oop_FNs[1:,7]]--- -- 
# axs4.hist(mcc_f1_ref_FNs_circ[1:,0], label='Reference Model', color='g') #, histtype='step', linewidth=2   
# axs4.legend()
# axs4.hist([mcc_f1_oop_FNs_circ[1:,1], mcc_f1_oop_FNs_circ[1:,2], mcc_f1_oop_FNs_circ[1:,3],mcc_f1_oop_FNs_circ[1:,4],
#             mcc_f1_oop_FNs_circ[1:,5], mcc_f1_oop_FNs_circ[1:,6]], 
#           label=['FN=2', 'FN=3', 'FN=4','FN=5','FN=6','FN=7'], alpha=0.3, color=['r', 'r', 'r','r','r','r'], 
#           linewidth=2, linestyle='--') #, roc_oop_FNs[1:,7]]--- -- 

# axs4.set_title('Circle Grid - MCC-F1 value: Reference Model vs Coulomb Model')
# axs4.set_xlabel('MCf1 value')


# AUC Circle Grid
nbins=50

fig5, axs5 = plt.subplots()
fig5.set_size_inches(8, 6)

axs5.plot(np.histogram(roc_oop_FNs_circ[1:,0], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_circ[1:,0], bins=nbins)[0]),
          color='r')  #,label='Coulomb Model' 
axs5.plot(np.histogram(roc_ref_FNs_circ[1:,0], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_circ[1:,0], bins=nbins)[0])),
            color='b') #, label='Reference Model'
axs5.legend(['Coulomb Model', 'Reference Model'], fontsize=14)

axs5.plot(np.histogram(roc_oop_FNs_circ[1:,1], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_circ[1:,1], bins=nbins)[0]),
          color='r',label='SE=1',alpha=0.8, linestyle='dotted') 
# axs5.plot(np.histogram(roc_ref_FNs_circ[1:,1], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_circ[1:,1], bins=nbins)[0])),
#           '--', color='g')


axs5.plot(np.histogram(roc_oop_FNs_circ[1:,2], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_circ[1:,1], bins=nbins)[0]),
          color='r',label='SE=2',alpha=0.8, linestyle='dotted') 
# axs5.plot(np.histogram(roc_ref_FNs_circ[1:,1], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_circ[1:,1], bins=nbins)[0])),
#           '--', color='g')

axs5.plot(np.histogram(roc_oop_FNs_circ[1:,3], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_circ[1:,2], bins=nbins)[0]), 
          color='r',label='SE=3', alpha=0.8, linestyle='dotted') 
# axs5.plot(np.histogram(roc_ref_FNs_circ[1:,2], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_circ[1:,2], bins=nbins)[0])),
#           '--', color='r')

axs5.plot(np.histogram(roc_oop_FNs_circ[1:,4], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_circ[1:,3], bins=nbins)[0]),
        color='r',label='SE=4', alpha=0.8, linestyle='dotted') 
# axs5.plot(np.histogram(roc_ref_FNs_circ[1:,3], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_circ[1:,3], bins=nbins)[0])),
#           '--', color='y')

axs5.plot(np.histogram(roc_oop_FNs_circ[1:,5], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_circ[1:,4], bins=nbins)[0]),
          color='r', label='SE=5', alpha=0.8, linestyle='dotted') 
# axs5.plot(np.histogram(roc_ref_FNs_circ[1:,4], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_circ[1:,4], bins=nbins)[0])),
#           '--', color='m')

axs5.plot(np.histogram(roc_oop_FNs_circ[1:,6], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_circ[1:,5], bins=nbins)[0]),
          color='r', label='SE=6', alpha=0.8, linestyle='dotted') 
# axs5.plot(np.histogram(roc_ref_FNs_circ[1:,5], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_circ[1:,5], bins=nbins)[0])),
#           '--', color='c')

axs5.plot(np.histogram(roc_oop_FNs_circ[1:,7], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_circ[1:,6], bins=nbins)[0]),
           color='r', label='SE=7', alpha=0.8, linestyle='dotted') 
# axs5.plot(np.histogram(roc_ref_FNs_circ[1:,6], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_circ[1:,6], bins=nbins)[0])),
#           '--', color='k')

labelLines(axs5.get_lines(), zorder=5, align = True)
#axs5.set_title('Circle Grid - AUC value: Reference Model vs Coulomb Model', fontsize=18)
axs5.set_xlabel('AUC values', fontsize=14)
axs5.set_ylabel('Distribution of of AUC values', fontsize=14)


# AUC Square grid
fig6, axs6 = plt.subplots()
fig6.set_size_inches(8, 6)

axs6.plot(np.histogram(roc_oop_FNs_square[1:,0], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_square[1:,0], bins=nbins)[0]),
          color='r')  #,label='Coulomb Model' 
axs6.plot(np.histogram(roc_ref_FNs_square[1:,0], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_square[1:,0], bins=nbins)[0])),
            color='b')
axs6.legend(['Coulomb Model', 'Reference Model'], fontsize=14)

axs6.plot(np.histogram(roc_oop_FNs_square[1:,1], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_square[1:,1], bins=nbins)[0]),
          color='r',label='SE=1', alpha=0.8, linestyle='dotted') 
# axs6.plot(np.histogram(roc_ref_FNs_square[1:,1], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_square[1:,1], bins=nbins)[0])),
#           '--', color='g')

axs6.plot(np.histogram(roc_oop_FNs_square[1:,2], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_square[1:,1], bins=nbins)[0]),
          color='r',label='SE=2', alpha=0.8, linestyle='dotted') 
# axs6.plot(np.histogram(roc_ref_FNs_square[1:,1], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_square[1:,1], bins=nbins)[0])),
#           '--', color='g')

axs6.plot(np.histogram(roc_oop_FNs_square[1:,3], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_square[1:,2], bins=nbins)[0]), 
          color='r',label='SE=3', alpha=0.8, linestyle='dotted') 
# axs6.plot(np.histogram(roc_ref_FNs_square[1:,2], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_square[1:,2], bins=nbins)[0])),
#           '--', color='r')

axs6.plot(np.histogram(roc_oop_FNs_square[1:,4], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_square[1:,3], bins=nbins)[0]),
        color='r',label='SE=4', alpha=0.8, linestyle='dotted') 
# axs6.plot(np.histogram(roc_ref_FNs_square[1:,3], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_square[1:,3], bins=nbins)[0])),
#           '--', color='y')

axs6.plot(np.histogram(roc_oop_FNs_square[1:,5], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_square[1:,4], bins=nbins)[0]),
          color='r', label='SE=5', alpha=0.8, linestyle='dotted') 
# axs6.plot(np.histogram(roc_ref_FNs_square[1:,4], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_square[1:,4], bins=nbins)[0])),
#           '--', color='m')

axs6.plot(np.histogram(roc_oop_FNs_square[1:,6], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_square[1:,5], bins=nbins)[0]),
          color='r', label='SE=6', alpha=0.8, linestyle='dotted') 
# axs6.plot(np.histogram(roc_ref_FNs_square[1:,5], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_square[1:,5], bins=nbins)[0])),
#           '--', color='c')

axs6.plot(np.histogram(roc_oop_FNs_square[1:,7], bins=nbins)[1][:-1], np.cumsum(np.histogram(roc_oop_FNs_square[1:,6], bins=nbins)[0]),
          color='r', label='SE=7', alpha=0.8, linestyle='dotted') 
# axs6.plot(np.histogram(roc_ref_FNs_square[1:,6], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(roc_ref_FNs_square[1:,6], bins=nbins)[0])),
#           '--', color='k')

labelLines(axs6.get_lines(), zorder=5, align = True)
#axs6.set_title('Square Grid - AUC value: Reference Model vs Coulomb Model', fontsize=18)
axs6.set_xlabel('AUC values', fontsize=14)
axs6.set_ylabel('Distribution of of AUC values', fontsize=14)


# Square Grid MCC-F1
fig7, axs7 = plt.subplots()
fig7.set_size_inches(8, 6)

axs7.plot(np.histogram(mcc_f1_oop_FNs_square[1:,0], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_square[1:,0], bins=nbins)[0]),
          color='r') #label='Coulomb Model' 
axs7.plot(np.histogram(mcc_f1_ref_FNs_square[1:,0], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_square[1:,0], bins=nbins)[0])),
          color='b') #label='Coulomb Model'
axs7.legend(['Coulomb Model', 'Reference Model'], fontsize=14)


axs7.plot(np.histogram(mcc_f1_oop_FNs_square[1:,1], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_square[1:,1], bins=nbins)[0]),
          color='r',label='SE=1', alpha=0.8, linestyle='dotted') 
# axs7.plot(np.histogram(mcc_f1_ref_FNs_square[1:,1], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_square[1:,1], bins=nbins)[0])),
#           '--', color='g')

axs7.plot(np.histogram(mcc_f1_oop_FNs_square[1:,2], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_square[1:,1], bins=nbins)[0]),
          color='r',label='SE=2', alpha=0.8, linestyle='dotted') 
# axs7.plot(np.histogram(mcc_f1_ref_FNs_square[1:,1], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_square[1:,1], bins=nbins)[0])),
#           '--', color='g')

axs7.plot(np.histogram(mcc_f1_oop_FNs_square[1:,3], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_square[1:,2], bins=nbins)[0]), 
          color='r',label='SE=3', alpha=0.8, linestyle='dotted') 
# axs7.plot(np.histogram(mcc_f1_ref_FNs_square[1:,2], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_square[1:,2], bins=nbins)[0])),
#           '--', color='r')

axs7.plot(np.histogram(mcc_f1_oop_FNs_square[1:,4], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_square[1:,3], bins=nbins)[0]),
        color='r',label='SE=4', alpha=0.8, linestyle='dotted') 
# axs7.plot(np.histogram(mcc_f1_ref_FNs_square[1:,3], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_square[1:,3], bins=nbins)[0])),
#           '--', color='y')

axs7.plot(np.histogram(mcc_f1_oop_FNs_square[1:,5], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_square[1:,4], bins=nbins)[0]),
          color='r', label='SE=5', alpha=0.8, linestyle='dotted') 
# axs7.plot(np.histogram(mcc_f1_ref_FNs_square[1:,4], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_square[1:,4], bins=nbins)[0])),
#           '--', color='m')

axs7.plot(np.histogram(mcc_f1_oop_FNs_square[1:,6], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_square[1:,5], bins=nbins)[0]),
          color='r', label='SE=6', alpha=0.8, linestyle='dotted') 
# axs7.plot(np.histogram(mcc_f1_ref_FNs_square[1:,5], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_square[1:,5], bins=nbins)[0])),
#           '--', color='c')

axs7.plot(np.histogram(mcc_f1_oop_FNs_square[1:,7], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_square[1:,6], bins=nbins)[0]),
          color='r', label='SE=7', alpha=0.8, linestyle='dotted') 
# axs7.plot(np.histogram(mcc_f1_ref_FNs_square[1:,6], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_square[1:,6], bins=nbins)[0])),
#           '--', color='k')

labelLines(axs7.get_lines(), zorder=5, align = True)
#axs7.set_title('Square Grid - MCC-F1 value: Reference Model vs Coulomb Model', fontsize=18)
axs7.set_xlabel('MCC-F1 metric', fontsize=14)
axs7.set_ylabel('Distribution of MCC-F1 metric', fontsize=14)


#Circle Grid --- MCC-F1
fig8, axs8 = plt.subplots()
fig8.set_size_inches(8, 6)

axs8.plot(np.histogram(mcc_f1_oop_FNs_circ[1:,0], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_circ[1:,0], bins=nbins)[0]),
          color='r')  #,label='Coulomb Model'
axs8.plot(np.histogram(mcc_f1_ref_FNs_circ[1:,0], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_circ[1:,0], bins=nbins)[0])),
            color='b') #label='Reference Model'
axs8.legend(['Coulomb Model', 'Reference Model'], fontsize=14)

axs8.plot(np.histogram(mcc_f1_oop_FNs_circ[1:,1], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_circ[1:,1], bins=nbins)[0]),
          color='r',label='SE=1', alpha=0.8, linestyle='dotted') 
# axs8.plot(np.histogram(mcc_f1_ref_FNs_circ[1:,1], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_circ[1:,1], bins=nbins)[0])),
#           '--', color='g')


axs8.plot(np.histogram(mcc_f1_oop_FNs_circ[1:,2], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_circ[1:,1], bins=nbins)[0]),
          color='r',label='SE=2', alpha=0.8, linestyle='dotted') 
# axs8.plot(np.histogram(mcc_f1_ref_FNs_circ[1:,1], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_circ[1:,1], bins=nbins)[0])),
#           '--', color='g')

axs8.plot(np.histogram(mcc_f1_oop_FNs_circ[1:,3], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_circ[1:,2], bins=nbins)[0]), 
          color='r',label='SE=3', alpha=0.8, linestyle='dotted') 
# axs8.plot(np.histogram(mcc_f1_ref_FNs_circ[1:,2], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_circ[1:,2], bins=nbins)[0])),
#           '--', color='r')

axs8.plot(np.histogram(mcc_f1_oop_FNs_circ[1:,4], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_circ[1:,3], bins=nbins)[0]),
        color='r',label='SE=4', alpha=0.8, linestyle='dotted') 
# axs8.plot(np.histogram(mcc_f1_ref_FNs_circ[1:,3], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_circ[1:,3], bins=nbins)[0])),
#           '--', color='y')

axs8.plot(np.histogram(mcc_f1_oop_FNs_circ[1:,5], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_circ[1:,4], bins=nbins)[0]),
          color='r', label='SE=5', alpha=0.8, linestyle='dotted') 
# axs8.plot(np.histogram(mcc_f1_ref_FNs_circ[1:,4], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_circ[1:,4], bins=nbins)[0])),
#           '--', color='m')

axs8.plot(np.histogram(mcc_f1_oop_FNs_circ[1:,6], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_circ[1:,5], bins=nbins)[0]),
          color='r', label='SE=6', alpha=0.8, linestyle='dotted') 
# axs8.plot(np.histogram(mcc_f1_ref_FNs_circ[1:,5], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_circ[1:,5], bins=nbins)[0])),
#           '--', color='c')

axs8.plot(np.histogram(mcc_f1_oop_FNs_circ[1:,7], bins=nbins)[1][:-1], np.cumsum(np.histogram(mcc_f1_oop_FNs_circ[1:,6], bins=nbins)[0]),
          color='r', label='SE=7', alpha=0.8, linestyle='dotted') 
# axs8.plot(np.histogram(mcc_f1_ref_FNs_circ[1:,6], bins=nbins)[1][:-1],np.flip(np.cumsum(np.histogram(mcc_f1_ref_FNs_circ[1:,6], bins=nbins)[0])),
#           '--', color='k')

labelLines(axs8.get_lines(), zorder=5, align = True)
#axs8.set_title('Circle Grid - MCC-F1 value: Reference Model vs Coulomb Model', fontsize=18)
axs8.set_xlabel('MCC-F1 metric', fontsize=14)
axs8.set_ylabel('Dstribution of MCC-F1 metric', fontsize=14)

fig5.tight_layout()
fig6.tight_layout()
fig7.tight_layout()
fig8.tight_layout()

if stress_MAS:
    stress_name = 'MAS'
else:
    stress_name = 'OOP'

if save_results:
    # axs1.figure.savefig('results/Square_grid_ROC_'+stress_name+'_vs_'+model_type+'_hist.svg')
    # axs2.figure.savefig('results/Circular_grid_ROC_'+stress_name+'_vs_'+model_type+'_hist.svg')
    # axs3.figure.savefig('results/Square_grid_MCC-F1_'+stress_name+'_vs_'+model_type+'_hist.svg')
    # axs4.figure.savefig('results/Circle_grid_MCC-F1_'+stress_name+'vs_'+model_type+'_hist.svg')
    
    axs5.figure.savefig('square_circle/Circular_grid_ROC_'+stress_name+'_:_'+model_type+'_cumulative.png', dpi=300)
    axs6.figure.savefig('square_circle/Square_grid_ROC_'+stress_name+'_:_'+model_type+'_cumulative.png', dpi = 300)
    axs7.figure.savefig('square_circle/Square_grid_MCC-F1_'+stress_name+'_:_'+model_type+'_cumulative.png', dpi=300)
    axs8.figure.savefig('square_circle/Circle_grid_MCC-F1_'+stress_name+'_:_'+model_type+'_cumulative.png', dpi=300)
    
    