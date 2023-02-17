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
# import matplotlib.pyplot as plt
# from shapely.geometry import box, shape, mapping, Point, Polygon
# import json
from grid_operations import bounds_to_polygons, create_square_grid_bounds, create_circular_grid, forecast_aggregation
from evaluations import calc_ROC_MCCF1, calc_mccf1_metric
from utils import write_geojson_feature, generate_gridded_eqs, plot_performance_distribution
import config

def main():
    print ('Generating Figure 2b, 2d, 3b, 3d .....')
    if config.aggregation:
        #---------------Start of code
        events = pandas.read_pickle('../data/calculated_stress.pkl')
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
        filename = '../data/model_grid_mas_rescaledCoords'
    #    write_geojson_feature(model_grid,stress_data[:,3],filename )
        ###--========-- Convert the forecast points into grid cells here .....
        if config.stress_MAS:
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
        fault_line = numpy.array([[0,yy] for yy in numpy.arange(-5,6,1)])
    
        dist_line = []
        # dist_point = []
        for coord in coords:
            dist_line.append(min(numpy.sqrt((coord[0]-fault_line[:,0])**2 + (coord[1]-fault_line[:,1])**2)))
            
        dist_line = numpy.array(dist_line)
        # dist_point = numpy.array(dist_point)
        c = 2 #Static stress constant
        static_stress = c* dist_line ** -2
       
        
        if numpy.isinf(max(static_stress)):
            static_stress[numpy.where(numpy.isinf(static_stress) == True)] = numpy.sort(static_stress)[-2]
    
        
        stress_data = numpy.column_stack((stress_data, static_stress))
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
    
    
    if config.stress_MAS:
            fname_square = '../data/square_grid_aggregated_stress_'+str(len(square_grid))+'_MAS.csv'
            fname_circle = '../data/circle_grid_aggregated_stress_'+str(len(circle_grid))+'_MAS.csv'
    else:
            fname_square = '../data/square_grid_aggregated_stress_'+str(len(square_grid))+'_OOP.csv'
            fname_circle = '../data/circle_grid_aggregated_stress_'+str(len(circle_grid))+'_OOP.csv'
    
    
    if config.aggregation:
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
    
    for i in range(config.num_sim):
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
                                                                                     model_type=config.model_for_sim, Num_eqs=[config.total_eqs-fn, fn])
            #write_geojson_feature(square_grid_arranged.tolist(), stress_data_square_arranged[:,2],
            #'results_new/square_grid_eqs_arranged_'+str(fn))
            #---
            circle_grid_arranged, stress_data_circle_arranged = generate_gridded_eqs(circle_grid, stress_data_circle, 
                                                                                     model_type=config.model_for_sim, Num_eqs=[config.total_eqs-fn, fn])
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
            
     
    # AUC Square grid
    
    axs1 = plot_performance_distribution(roc_oop_FNs_square, roc_ref_FNs_square)
    axs1.set_xlabel('AUC values', fontsize=14)
    axs1.set_ylabel('Distribution of of AUC values', fontsize=14)
    
    # AUC Circle Grid
    # nbins=50
    axs2 = plot_performance_distribution(roc_oop_FNs_circ, roc_ref_FNs_circ)
    axs2.set_xlabel('AUC values', fontsize=14)
    axs2.set_ylabel('Distribution of of AUC values', fontsize=14)
    
    # Square Grid MCC-F1
    axs3 = plot_performance_distribution(mcc_f1_oop_FNs_square, mcc_f1_ref_FNs_square)
    axs3.set_xlabel('MCC-F1 metric', fontsize=14)
    axs3.set_ylabel('Distribution of MCC-F1 metric', fontsize=14)
    
    #Circle Grid --- MCC-F1
    axs4 = plot_performance_distribution(mcc_f1_oop_FNs_circ, mcc_f1_ref_FNs_circ)
    axs4.set_xlabel('MCC-F1 metric', fontsize=14)
    axs4.set_ylabel('Dstribution of MCC-F1 metric', fontsize=14)
    
    axs1.figure.tight_layout()
    axs2.figure.tight_layout()
    axs3.figure.tight_layout()
    axs4.figure.tight_layout()
    
    if config.stress_MAS:
        stress_name = 'MAS'
    else:
        stress_name = 'OOP'
    
    if config.save_results:
        axs1.figure.savefig('../output/Fig2b_Square_grid_ROC_'+stress_name+'_:_'+config.model_for_sim+'_cumulative.png', dpi = 300) #fig2b
        axs2.figure.savefig('../output/Fig3b_Circular_grid_ROC_'+stress_name+'_:_'+config.model_for_sim+'_cumulative.png', dpi=300) #fig3b
        axs3.figure.savefig('../output/Fig2d_Square_grid_MCC-F1_'+stress_name+'_:_'+config.model_for_sim+'_cumulative.png', dpi=300) #fig2d
        axs4.figure.savefig('../output/Fig3d_Circle_grid_MCC-F1_'+stress_name+'_:_'+config.model_for_sim+'_cumulative.png', dpi=300) #fig3d
    

if __name__ == "__main__":
    
    main()