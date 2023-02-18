#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 17:05:38 2022

@author: khawaja
"""
import numpy
import pandas
#from utils import grid_origin_coordinates, grid_top_right_coordinates #, geographical_area_from_bounds, write_geojson_feature
from experiment.grid_operations import locate_eq_to_grid3D #forecast_aggregation_bounds, bounds_to_polygons, 
from experiment.evaluations import calc_ROC_MCCF1, calc_mccf1_metric
import matplotlib.pyplot as plt
import experiment.config as config
#from shapely.geometry import Point

# which_earthquake =  'chichi'  #'landers' # #  #'

grids = config.grids
def evaluate_real_aftershock_model(earthquake_name):
    
    if earthquake_name == 'chichi':
        eq_name ='s1999CHICHI01MAxx' 
        loc = [120.840,23.869, 7.0]
        Mc = 3.0
    elif earthquake_name == 'landers':
        eq_name = 's1992LANDER01WALD'
        # folder = 'landers/'
        loc = [-116.430, 34.200, 7.0]
        Mc = 2.0
    else:
        print('select the right earthquake')
    
    # model = 'OOP' #'MAS' # #
    auc_cfs = []
    auc_ref = []
    mcc_f1_cfs = []
    mcc_f1_ref = []
    folder ='data/' 

    
    display_depth=7.5
    #Get the catalog and filter for depth.
    cat = numpy.loadtxt(folder+'Aftershocks_'+eq_name+'365.csv', skiprows=1, delimiter=',') #  
    #To ensure cut-off mag
    cat = cat[cat[:,2]>=Mc]
    #
    # cat = cat[cat[:,3]<=8]
    # cat = cat[cat[:,3]>=6]
    # numpy.savetxt(folder+folder[:-1]+'_cat_L14km.csv', cat, delimiter=',', header='lat,lon,mag,depth', comments='')
    
    # colors=['black','red','blue','green']
    #1- Get the stress data (Chihi, Kashmir, or anyother)
    # ------stress for a grid...
    
    # fig1, axs1 = plt.subplots()
    # fig1.set_size_inches(8, 6)
    # axs1.set_title('ROC:Quadtree')
    # fig2, axs2 = plt.subplots()
    # fig2.set_size_inches(8, 6)
    # axs2.set_title('MCC-F1: Quadtree')
    
    for i, grid in enumerate(grids):
        stress_data = pandas.read_pickle(folder+eq_name+grid+'.pkl') 
        stress = stress_data.to_numpy()
        
        ##-------For LANDERS Depth Filter.....
        if earthquake_name == 'landers':
            stress = stress[stress[:,5] <=32]
        
        print(grid)
        # print(colors[i])
        # print('Stress len :',len(stress))
    
    #2- Get a grid
    
        qk = numpy.loadtxt(folder+eq_name+grid+'.txt', dtype='str')
        print('Quadkey len :', len(qk))
        # grid_bounds = numpy.column_stack((grid_origin_coordinates(qk), grid_top_right_coordinates(qk)))
        # grid_poly = bounds_to_polygons(grid_bounds)
         
        #Rescaled by volume to display .....
        # ---cell_volume = [quadtree_FuncFile.geographical_area_from_bounds(bb[0],bb[1],bb[2],bb[3]) * (bb[5]-bb[4]) for bb in grid_3d]
        
        #3- Aggregate them onto suitable grids (square, circular, quadtree) and evaluate. 
        
        # stress_loc = stress[:,:2]
        if config.stress_MAS:
            stress_data = stress[:,6]
        else:
            stress_data = stress[:,7]
            
        # ---its already aggregated now. So we dont need it anymore. 
        # # stress_aggregate = forecast_aggregation_bounds(grid_bounds, stress_loc, stress_data)
        #Plot for a particular depth
        # selected_depth = numpy.logical_and(stress[:,4]<=display_depth, display_depth<stress[:,5])
        # selected_stress = stress_data[selected_depth]
        
        # write_geojson_feature(grid_poly,selected_stress, folder+eq_name+'_'+model+'_'+str(display_depth)+grid) 
        
        #Save stress per unit volume (normaize er volume).
        # selected_data = stress[selected_depth]
        # cell_volume = [geographical_area_from_bounds(bb[0],bb[1],bb[2],bb[3]) * (bb[5]-bb[4]) for bb in selected_data]
        # stress_volume = selected_stress/cell_volume
        # write_geojson_feature(grid_poly,selected_stress/cell_volume, folder+eq_name+'_'+model+'_'+str(display_depth)+grid+'_per_volume') 
        
        
        #3a Get the REF forecast
        # ---Lon, Lat, Depth 
        center_coords = numpy.column_stack(( numpy.column_stack(((stress[:,0]+stress[:,2])/2, 
                                         (stress[:,1]+stress[:,3])/2)), (stress[:,4]+stress[:,5])/2))
        
        dist = []
        for coord in center_coords:
    #        dist.append(numpy.sqrt((coord[0]-chichi_loc[0])**2 + (coord[1]-chichi_loc[1])**2) + (coord[2]-chichi_loc[2])**2 )
            dist.append(numpy.sqrt((coord[0]-loc[0])**2 + (coord[1]-loc[1])**2) + (coord[2]-loc[2])**2 )
                
        dist = numpy.array(dist) 
        c = 2 #Static stress constant
        static_stress = c* dist ** -2
        #static_grid_stress = numpy.column_stack((grid_bounds, static_stress))
        #    
        #if numpy.isinf(max(static_stress)):
        #    static_grid_stress = numpy.delete(static_grid_stress, numpy.where(numpy.isinf(static_stress) == True),0)
        # ---For a particular depth.
        # write_geojson_feature(grid_poly,static_stress[selected_depth], folder+eq_name+'_REF'+'_'+str(display_depth)+grid) 
        
        #4- Get the aftershocks
    
        # stress_coords = stress[:,:6]    
        gridded_cat = locate_eq_to_grid3D(stress[:,:6], cat)
        print('------Total Gridded EQs :', sum(gridded_cat))
        
        print('-----The Percentage of active cells :', len(gridded_cat[gridded_cat>0]) / len(gridded_cat) *100)
    
        #5- Run the test
        oop_ROCdata, oop_MCC_F1 = calc_ROC_MCCF1(stress_data, gridded_cat)
        fpr = oop_ROCdata[:,0]
        tpr = oop_ROCdata[:,1]
        auc = numpy.round(abs(numpy.trapz(tpr, fpr)), 3)
        f1 = oop_MCC_F1[:,0]
        mcc = oop_MCC_F1[:,1]
        mcc_f1 = numpy.round(calc_mccf1_metric(f1,mcc, min_dist=True),3)
        
        auc_cfs.append(auc)
        mcc_f1_cfs.append(mcc_f1)
        
        
        # axs1.plot(fpr, tpr, label = model+grid+': '+str(auc), color=colors[i])
        # axs2.plot(f1,mcc, label = model+grid+': '+str(mcc_f1), color=colors[i])
        
        #Ref 
        oop_ROCdata, oop_MCC_F1 = calc_ROC_MCCF1(static_stress, gridded_cat)
        fpr = oop_ROCdata[:,0]
        tpr = oop_ROCdata[:,1]
        auc = numpy.round(abs(numpy.trapz(tpr, fpr)), 3)
        f1 = oop_MCC_F1[:,0]
        mcc = oop_MCC_F1[:,1]
        mcc_f1 = numpy.round(calc_mccf1_metric(f1,mcc, min_dist=True),3)
        
        auc_ref.append(auc)
        mcc_f1_ref.append(mcc_f1)

    
    grid_label = ['L12', 'L13', 'L14','Multi_res'] #' 
    fig3, axs3 = plt.subplots()
    fig3.set_size_inches(8, 6)
    
    axs3.plot(grid_label, auc_cfs, color ='red', label='$\Delta CFS$ : AUC')
    axs3.plot(grid_label, auc_ref, color='blue', label='R : AUC')
    
    axs3.plot(grid_label, mcc_f1_cfs, '--', color ='red', label='$\Delta CFS$ : MCC-F$_1$')
    axs3.plot(grid_label, mcc_f1_ref, '--', color ='blue', label='R : MCC-F$_1$')
    
    axs3.legend(fontsize=14)
    axs3.figure.tight_layout()
    axs3.figure.savefig('output/'+eq_name+'_'+model+'_performance.png', dpi=300)
    # axs3.figure.savefig(folder+eq_name+'_'+model+'_performance.svg')


def main():
    
    
    earthquakes = ['chichi', 'landers']
    
    for earthquake_name in earthquakes:
        print('Earthquake :', earthquake_name)
        evaluate_real_aftershock_model(earthquake_name)
        
if __name__ == "__main__": 
    main()