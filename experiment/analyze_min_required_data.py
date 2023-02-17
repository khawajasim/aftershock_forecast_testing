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
from experiment.grid_operations import forecast_aggregation_mean, create_square_grid_bounds #, create_circular_grid, bounds_to_polygons
from experiment.evaluations import calc_ROC_MCCF1, calc_mccf1_metric
from experiment.utils import generate_active_cells #generate_gridded_eqs,  write_geojson_feature,
import experiment.config as config
# from labellines import labelLines


def main():
    print('Generating Figure 4')

    depth = 7 #5 if Event_3
    # num_sim = 10
    
    total_eqs = numpy.arange(15,321,2) #2485  # 
    
    
    
    if config.aggregation:
        #---------------Start of code
        events = pandas.read_pickle('data/calculated_stress.pkl')
        stress_data = events.to_numpy()
        
        stress_data = stress_data[stress_data[:,2] == depth]
        
        #Rescaling the degrees to the km distance --- -150 to 150
        stress_data[:,0] = numpy.round(stress_data[:,0]/max(stress_data[:,0]) * 150)
        stress_data[:,1] = numpy.round(stress_data[:,1]/max(stress_data[:,1]) * 150)
        
        #Filter within 100km x 100km
        within_range = numpy.logical_and(numpy.logical_and(stress_data[:,0]>=-100, stress_data[:,1]>=-100),
                          numpy.logical_and(stress_data[:,0]<100, stress_data[:,1]<100))
        
        stress_data= stress_data[within_range]
        
        #####-========------ Convert the forecast points into grid cells: MAKE A GEOJSON for QGIS----
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
        dist = []
        for coord in coords:
            dist.append(numpy.sqrt(coord[0]**2 + coord[1]**2))
            
        dist = numpy.array(dist) 
        c = 2 #Static stress constant
        static_stress = c* dist ** -2  #-3 (2 by Hardebeck)
        
        if numpy.isinf(max(static_stress)):
            static_stress[numpy.where(numpy.isinf(static_stress) == True)] = numpy.sort(static_stress)[-2]
    #        stress_data = numpy.delete(stress_data, numpy.where(numpy.isinf(static_stress) == True),0)
    #        cell_bounds = numpy.delete(cell_bounds, numpy.where(numpy.isinf(static_stress) == True),0)
        
        stress_data = numpy.column_stack((stress_data, static_stress))   
    
    #-----SQUARE TEST GRID
    square_bounds = create_square_grid_bounds([-100,-100], [100,100],40)
    
    # square_grid_aggregated_stress_1600_OOP
    if config.stress_MAS:
            fname_square = 'data/square_grid_aggregated_stress_'+str(len(square_bounds))+'_MAS.csv'
    #        fname_circle = 'forecast_data/circle_grid_aggregated_stress_'+str(len(circle_grid))+'_MAS.csv'
    else:
            fname_square = 'data/square_grid_aggregated_stress_'+str(len(square_bounds))+'_OOP.csv'
    #        fname_circle = 'forecast_data/circle_grid_aggregated_stress_'+str(len(circle_grid))+'_OOP.csv'
    
    
    if config.aggregation:
    #    square_grid_oop = forecast_aggregation(square_grid, model_grid, stress_data[:,3])
        square_grid_cfs = forecast_aggregation_mean(square_bounds, stress_data[:,:2], stress_data[:,3])
        # write_geojson_feature(square_grid,square_grid_cfs,'forecast_data/square_grid_mas_aggregated_meanArea')
        
    #    square_grid_ref = forecast_aggregation(square_grid, model_grid, stress_data[:,4])
        square_grid_ref = forecast_aggregation_mean(square_bounds, stress_data[:,:2], stress_data[:,4])
        # write_geojson_feature(square_grid, square_grid_ref, 'forecast_data/square_grid_ref_aggregated_meanArea')
        
        stress_data_square = numpy.column_stack((square_grid_cfs, square_grid_ref))
        
    
        numpy.savetxt(fname_square, 
                      stress_data_square, delimiter=',', header='CFS,Ref', comments='')
    #    numpy.savetxt(fname_circle, 
    #                  stress_data_circle, delimiter=',', header='MAS,Ref', comments='')
    
    #------ Locate earthquakes to grid ---
    # numpy.savetxt('stress_data_square_test_grid_'+str(len(stress_data_square))+'.csv', stress_data_square, delimiter=',', header='oop,ref', comments='')
    # numpy.savetxt('stress_data_circle_test_grid_'+str(len(stress_data_circle))+'.csv', stress_data_circle, delimiter=',', header='oop,ref', comments='')
    
    stress_data_square = numpy.loadtxt(fname_square , delimiter=',', skiprows = 1)
    
    #Generate a stress_data variable for use in function "generate_active_cells"
            #
            #stress_data: n x 5: [x, y, depth, coulomb_stress, reference_stress]
    stress = numpy.column_stack((numpy.column_stack((square_bounds[:,:2], numpy.ones(len(square_bounds))*depth)),
                                stress_data_square))
    
    
    
    #FNs = [0,1,2,3,4,5,6,7]
    FNs = [0,10]
    data_roc = {}
    data_mccf1 = {}
    
    for fn in FNs:
        print('Number of False Negative :',fn)
        roc_percentile_diff = []
        mccf1_percentile_diff = []
        for eqs in total_eqs:   
            roc_cfs_square = []
            roc_ref_square = []
            mcc_f1_cfs_square = []
            mcc_f1_ref_square = []
            
            for i in range(config.num_sim):
        
        #        for fn in FNs:
                
                stress_eqs = generate_active_cells(stress, Num_eqs=[eqs-fn, fn],  model_type= config.model_for_sim) #
            
                #----============== PLOTS for Square Grids.
                oop_stress = stress_eqs[:,3]
                static_stress = stress_eqs[:,4]
                earthquakes = stress_eqs[:,5]
                # print('No. of Active cells for SQUARE Grid:', len(earthquakes[earthquakes>0]))
                # print('----Active FNs cells in SQUARE Grid:', len(oop_stress[numpy.logical_and(earthquakes>0, oop_stress<0)]))
                
                #--- ------Plot ROC OOP stress
                oop_ROCdata, oop_MCC_F1 = calc_ROC_MCCF1(oop_stress, earthquakes)
                fpr = oop_ROCdata[:,0]
                tpr = oop_ROCdata[:,1]
                auc = numpy.round(abs(numpy.trapz(tpr, fpr)), 3)
                roc_cfs_square.append(auc)
                    
                f1 = oop_MCC_F1[:,0]
                mcc = oop_MCC_F1[:,1]
                mcc_f1 = numpy.round(calc_mccf1_metric(f1,mcc, min_dist=True),3)
                mcc_f1_cfs_square.append(mcc_f1)
                
                static_ROCdata, static_MCC_F1 = calc_ROC_MCCF1(static_stress, earthquakes)
                fpr = static_ROCdata[:,0]
                tpr = static_ROCdata[:,1]
                auc = numpy.round(abs(numpy.trapz(tpr, fpr)), 3)
                roc_ref_square.append(auc)
                
                f1 = static_MCC_F1[:,0]
                mcc = static_MCC_F1[:,1]
                mcc_f1 = calc_mccf1_metric(f1,mcc, min_dist=True)
                mcc_f1_ref_square.append(mcc_f1)
        
            percentile_diff = numpy.percentile(roc_cfs_square, 1) - numpy.percentile(roc_ref_square, 99)
            roc_percentile_diff.append([eqs, percentile_diff])
            
            percentile_diff = numpy.percentile(mcc_f1_cfs_square,1) - numpy.percentile(mcc_f1_ref_square,99)
            mccf1_percentile_diff.append([eqs, percentile_diff])
        
        data_roc[fn] = roc_percentile_diff
        data_mccf1[fn] = mccf1_percentile_diff

    
    fig1, axs1 = plt.subplots()
    fig1.set_size_inches(8, 6)
    
    fn=0
    roc_perf = numpy.array(data_roc[fn])
    mccf1_perf = numpy.array(data_mccf1[fn])
    
    #roc_percentile_diff = numpy.array(roc_percentile_diff)
    axs1.plot(roc_perf[:,0]/len(stress)*100, roc_perf[:,1], label='ROC (SE=0)', color='black')
    #mccf1_percentile_diff = numpy.array(mccf1_percentile_diff)
    axs1.plot(mccf1_perf[:,0]/len(stress)*100, mccf1_perf[:,1], label='MCC-F1 (SE=0)', color='blue')
    
    fn=10
    roc_perf = numpy.array(data_roc[fn])
    mccf1_perf = numpy.array(data_mccf1[fn])
    
    axs1.plot(roc_perf[:,0]/len(stress)*100, roc_perf[:,1], label='ROC (SE=10)', color='black', linestyle='dotted')
    #mccf1_percentile_diff = numpy.array(mccf1_percentile_diff)
    axs1.plot(mccf1_perf[:,0]/len(stress)*100, mccf1_perf[:,1], label='MCC-F1 (SE=10)', color='blue', linestyle='dotted')
    axs1.axhline(y=0, color ='red', linestyle='--')
    
    
    
    axs1.legend(fontsize=14)
    axs1.set_xlabel('Percentage of Active cells (%)', fontsize=14)
    # axs1.set_ylabel('Perforamnce difference: $\delta$CFS($1^{th}$) - R($99^{th}$)', fontsize=14)
    axs1.set_ylabel('Quantile distance: $\Delta$ CFS(1\%) - R(99\%)', fontsize=14)
    fig1.tight_layout()
    axs1.figure.savefig('output/active_cells_ratio_FN=0_10.png', dpi=300)
    axs1.figure.savefig('output/active_cells_ratio_FN=0_10.svg')


if __name__ == "__main()__":
    main()

    