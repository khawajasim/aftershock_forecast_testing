#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 16:08:12 2022

@author: khawajasim
"""
# import scipy
import numpy
# import numpy as np
import pandas
import matplotlib.pyplot as plt
# from shapely.geometry import box, shape, mapping, Point, Polygon
# import json
from grid_operations import bounds_to_polygons, create_square_grid_bounds, create_circular_grid, forecast_aggregation
from evaluations import calc_ROC_MCCF1, calc_mccf1_metric
from utils import write_geojson_feature, generate_gridded_eqs, plot_grid_polygon



introduce_FN = True
aggregation = False
save_results = False
stress_MAS= True   #True for MAS and False for OOP      
model_type =  'coulomb' #'ref' 
total_eqs = 400 #2485  # #

#---------------Start of code
if aggregation:
    #Running code from scratch and loading high resolution forecast and then aggregate it on circle and square. 
    events = pandas.read_pickle('coulomb_forecast/Event_4.pkl')
    stress_data = events.to_numpy()
    depth = 7 #5 if Event_3
    stress_data = stress_data[stress_data[:,2] == depth]
    
    #Rescaling the (x,y) coordinates, just for better readibility. -150 to 150
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
    
    ####--------NEXT STEPS FROM HERE------ (10-10-2022)
    # 1 -  Generate distance based stress
    
    #Choose location data and OOP stress. #4=OOP, 3=MAS
    if stress_MAS:
        stress_data = numpy.column_stack((stress_data[:,:3], stress_data[:,3])) 
    else:
        stress_data = numpy.column_stack((stress_data[:,:3], stress_data[:,4]))
    
    # # Deleted the centeral point of the whole grid, to make it symetric for plotting with "plot_stress" function.
    # # And secondly to calculate the static stress from origin point (0,0). Because distance from origin point is zero, leading to inf stress.
    # stress_data = numpy.delete(stress_data, numpy.where(numpy.logical_or(stress_data[:,0]==0, stress_data[:,1]==0) == True), 0 ) 
    
    #data = numpy.loadtxt('stress_oop.csv', skiprows=1, delimiter=',')
    
     #-- Generate Static Stress
     # center = [0,0]
    coords = stress_data[:,:2]
    fault_line = numpy.array([[0,yy] for yy in numpy.arange(-5,6,1)])

    dist_line = []
    # dist_point = []
    for coord in coords:
        dist_line.append(min(numpy.sqrt((coord[0]-fault_line[:,0])**2 + (coord[1]-fault_line[:,1])**2)))
        
    dist_line = numpy.array(dist_line)
    # dist_point = numpy.array(dist_point)
    c = 2 #Static stress constant
    static_stress = c* dist_line ** -2  #-3 (2 by Hardebeck)
    # stress_data = numpy.column_stack((stress_data, static_stress))
    
    if numpy.isinf(max(static_stress)):
        static_stress[numpy.where(numpy.isinf(static_stress) == True)] = numpy.sort(numpy.unique(static_stress))[-2]
        # stress_data = numpy.delete(stress_data, numpy.where(numpy.isinf(static_stress) == True),0)
        # cell_bounds = numpy.delete(cell_bounds, numpy.where(numpy.isinf(static_stress) == True),0)
        
    
    stress_data = numpy.column_stack((stress_data, static_stress))
    model_grid = bounds_to_polygons(cell_bounds)
    filename = 'square_circle/model_grid_ref'
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

ax = plot_grid_polygon(circle_grid)
ax.figure.savefig('square_circle/radial_grid.png', dpi=300)

if stress_MAS:
        fname_square = 'forecast_data/square_grid_aggregated_stress_'+str(len(square_grid))+'_MAS.csv'
        fname_circle = 'forecast_data/circle_grid_aggregated_stress_'+str(len(circle_grid))+'_MAS.csv'
else:
        fname_square = 'forecast_data/square_grid_aggregated_stress_'+str(len(square_grid))+'_OOP.csv'
        fname_circle = 'forecast_data/circle_grid_aggregated_stress_'+str(len(circle_grid))+'_OOP.csv'

if aggregation:
    # ---------Aggregate on Square grid
    square_grid_oop = forecast_aggregation(square_grid, model_grid, stress_data[:,3])
    # write_geojson_feature(square_grid,square_grid_oop,'square_circle/square_grid_mas_'+str(numpy.size(square_grid)))
    square_grid_ref = forecast_aggregation(square_grid, model_grid, stress_data[:,4])
    # write_geojson_feature(square_grid, square_grid_ref, 'square_circle/square_grid_ref_'+str(numpy.size(square_grid)))
    stress_data_square = numpy.column_stack((square_grid_oop, square_grid_ref))
    
    # -------Aggregate on circular grid
    circle_grid_oop = forecast_aggregation(circle_grid, model_grid, stress_data[:,3])
    # write_geojson_feature(circle_grid, circle_grid_oop,'square_circle/circle_grid_mas_'+str(numpy.size(circle_grid)))
    circle_grid_ref = forecast_aggregation(circle_grid, model_grid, stress_data[:,4])
    # write_geojson_feature(circle_grid, circle_grid_ref,'square_circle/circle_grid_ref_'+str(numpy.size(circle_grid)))
    stress_data_circle = numpy.column_stack((circle_grid_oop, circle_grid_ref))
    
    #------ Locate earthquakes to grid ---
    
    
    numpy.savetxt(fname_square, 
                  stress_data_square, delimiter=',', header='Stress,Ref', comments='')
    numpy.savetxt(fname_circle, 
                  stress_data_circle, delimiter=',', header='Stress,Ref', comments='')


stress_data_square = numpy.loadtxt(fname_square , delimiter=',', skiprows = 1)
stress_data_circle = numpy.loadtxt(fname_circle, delimiter=',', skiprows =1)



#-- Generate synthetic catalog
#Use probability map of static stress to generate catalog. 
#But generate events in only those cells which have positive OOP.

#The number obtained by searching http://www.isc.ac.uk/iscbulletin/search/catalogue/
#for 2 months after Amaratice earthquake


#--- Plot ROC OOP stress
fig1, axs1 = plt.subplots()
fig1.set_size_inches(8, 6)
#axs1.set_title('ROC: Square Grid')
axs1.set_xlabel('False Positive Rate', fontsize=14)
axs1.set_ylabel('True Positive Rate', fontsize=14)

fig2, axs2 = plt.subplots()
fig2.set_size_inches(8, 6)
#axs2.set_title('MCC-F1: Square Grid')
axs2.set_xlabel('F1', fontsize=14)
axs2.set_ylabel('MCC', fontsize=14)

fig3, axs3 = plt.subplots()
fig3.set_size_inches(8, 6)
#axs3.set_title('ROC: Circular Grid')
axs3.set_xlabel('False Positive Rate', fontsize=14)
axs3.set_ylabel('True Positive Rate', fontsize=14)

fig4, axs4 = plt.subplots()
fig4.set_size_inches(8, 6)
#axs4.set_title('MCC-F1: Circular Grid')
axs4.set_xlabel('False Positive Rate', fontsize=14)
axs4.set_ylabel('True Positive Rate', fontsize=14)


#for grid in grid_types:
    #---- FNs = 6
FNs = [0,1, 2,3,4,5,6,7]
for fn in FNs:
        print('False Negative :',fn)

        square_grid_arranged, stress_data_square_arranged = generate_gridded_eqs(square_grid, stress_data_square,
                                                                                model_type=model_type, Num_eqs=[total_eqs-fn, fn])
#        write_geojson_feature(square_grid_arranged.tolist(), stress_data_square_arranged[:,2],
#                              'square_circle/square_grid_eqs_arranged_'+str(fn)+'_'+model_type)
#        #---
        circle_grid_arranged, stress_data_circle_arranged = generate_gridded_eqs(circle_grid, stress_data_circle, 
                                                                                 model_type=model_type, Num_eqs=[total_eqs-fn, fn])
#        write_geojson_feature(circle_grid_arranged.tolist(), stress_data_circle_arranged[:,2],
#                              'square_circle/circle_grid_eqs_arranged_'+str(fn)+'_'+model_type)

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
        if fn == 0:
            axs1.plot(fpr, tpr, label = 'Coulomb Model AUC='+str(auc), color='red',alpha=1)
        elif fn == FNs[-1]:
            axs1.plot(fpr, tpr, color='red', alpha=0.20) #label = 'Coulomb Model AUC='+str(auc), 
        else:
            axs1.plot(fpr, tpr,  color='red',alpha=0.20) #label = 'OOP AUC FN_'+str(fn)+'='+str(auc),
        
        f1 = oop_MCC_F1[:,0]
        mcc = oop_MCC_F1[:,1]
        mcc_f1 = numpy.round(calc_mccf1_metric(f1,mcc, min_dist=True),3)

        if fn ==0:
            axs2.plot(f1,mcc, label = 'Coulomb Model MCC-F1='+str(mcc_f1), color='red', alpha=1)
        elif fn == FNs[-1]:
            axs2.plot(f1, mcc, color='red', alpha=0.20) #label = 'Coulomb Model MCC-F1='+str(mcc_f1),
        else:
            axs2.plot(f1,mcc,  color='red', alpha=0.20) #label = 'OOP MCC-F1 FN_'+str(fn)+'='+str(mcc_f1),
        
        
        # ---------Plot ROC for Reference Stress
        static_ROCdata, static_MCC_F1 = calc_ROC_MCCF1(static_stress, earthquakes)
        fpr = static_ROCdata[:,0]
        tpr = static_ROCdata[:,1]
        auc = numpy.round(abs(numpy.trapz(tpr, fpr)), 3)
        if fn ==0:
            axs1.plot(fpr, tpr, label = 'Reference Model AUC='+str(auc), color='blue', alpha=1)
        # elif fn == FNs[-1]:
        #     axs1.plot(fpr, tpr, label = 'Reference Model AUC='+str(auc), color='blue', alpha=0.20)
        # else:
        #     axs1.plot(fpr, tpr,color='blue', alpha=0.20) # label = 'Reference AUC FN_'+str(fn)+'='+str(auc), 
            
        f1 = static_MCC_F1[:,0]
        mcc = static_MCC_F1[:,1]
        mcc_f1 = numpy.round(calc_mccf1_metric(f1,mcc, min_dist=True),3)
        
        if fn == 0:
            axs2.plot(f1, mcc, label = 'Reference Model MCC-F1='+str(mcc_f1), color='blue', alpha=1)
        # elif fn == FNs[-1]:
        #     axs2.plot(f1, mcc, label = 'Reference Model MCC-F1='+str(mcc_f1), color='blue', alpha=0.20)
        # else:
        #     axs2.plot(f1, mcc, color='blue', alpha=0.20) #label = 'Reference MCC-F1 FN_'+str(fn)+'='+str(mcc_f1), 

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
        if fn == 0:
            axs3.plot(fpr, tpr, label = 'Coulomb Model AUC = '+str(auc), color='red',alpha=1)
        elif fn == FNs[-1]:
            axs3.plot(fpr, tpr, color='red', alpha=0.20) #label = 'Coulomb Model AUC = '+str(auc), 
        else:
            axs3.plot(fpr, tpr, color='red',alpha=0.20) #label = 'Coulomb AUC FN_'+str(fn)+'='+str(auc), 
        # if fn == max(FNs):
        #     axs1.plot(fpr, tpr, label = 'OOP AUC FN_'+str(fn)+'='+str(auc), color='black')

            
        
        f1 = oop_MCC_F1[:,0]
        mcc = oop_MCC_F1[:,1]
        mcc_f1 = numpy.round(calc_mccf1_metric(f1,mcc, min_dist=True),3)

        if fn ==0:
            axs4.plot(f1,mcc, label = 'Coulomb Model MCC-F1 = '+str(mcc_f1), color='red', alpha=1)
        elif fn == FNs[-1]:
            axs4.plot(f1, mcc,  color='red', alpha=0.20) #label = 'Coulomb Model MCC-F1 = '+str(mcc_f1),
        else:
            axs4.plot(f1,mcc,  color='red', alpha=0.20) #label = 'Coulomb MCC-F1 FN_'+str(fn)+'='+str(mcc_f1),

        
        # ---------Plot ROC for Reference Stress
        static_ROCdata, static_MCC_F1 = calc_ROC_MCCF1(static_stress, earthquakes)
        fpr = static_ROCdata[:,0]
        tpr = static_ROCdata[:,1]
        auc = numpy.round(abs(numpy.trapz(tpr, fpr)), 3)
        if fn ==0:
            axs3.plot(fpr, tpr, label = 'Reference Model AUC = '+str(auc), color='blue', alpha=1)
        # elif fn == FNs[-1]:
        #     axs3.plot(fpr, tpr, label = 'Reference Model AUC = '+str(auc), color='blue', alpha=0.20)
        # else:
        #     axs3.plot(fpr, tpr,  color='blue', alpha=0.20) #label = 'Reference AUC FN_'+str(fn)+'='+str(auc),
            
        f1 = static_MCC_F1[:,0]
        mcc = static_MCC_F1[:,1]
        mcc_f1 = numpy.round(calc_mccf1_metric(f1,mcc, min_dist=True),3)
        
        if fn == 0:
            axs4.plot(f1, mcc, label = 'Reference Model MCC-F1 = '+str(mcc_f1), color='blue', alpha=1)
        # elif fn == FNs[-1]:
        #     axs4.plot(f1, mcc, label = 'Reference Model MCC-F1 = '+str(mcc_f1), color='blue', alpha=0.20)
        # else:
        #     axs4.plot(f1, mcc,  color='blue', alpha=0.20) #label = 'Reference MCC-F1 FN_'+str(fn)+'='+str(mcc_f1),

        
axs1.legend(fontsize=14)
axs2.legend(fontsize=14)
axs3.legend(fontsize=14)
axs4.legend(fontsize=14)

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()


if stress_MAS:
    stress_name = 'MAS'
else:
    stress_name = 'OOP'
    
if save_results:
    axs1.figure.savefig('square_circle/Square_grid_ROC_'+stress_name+'_:_'+model_type+'.png', dpi=300) #Fig2a
    axs2.figure.savefig('square_circle/Square_grid_MCC-F1_'+stress_name+'_:_'+model_type+'.png',dpi=300) #Fig2c
    axs3.figure.savefig('square_circle/Circular_grid_ROC_'+stress_name+'_:_'+model_type+'.png', dpi=300) #Fig3a
    axs4.figure.savefig('square_circle/Circular_grid_MCC-F1_'+stress_name+'_:_'+model_type+'.png', dpi=300) #Fig3c