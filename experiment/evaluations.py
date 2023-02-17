#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 17:19:35 2022

@author: khawajasim
"""
import numpy


def calc_ROC_MCCF1(stress, earthquakes):
    ROCdata = []
    MCC_F1 = []
    for CFSth in numpy.sort(stress):
#        print(CFSth)
        TP = len(stress[((stress>=CFSth) & (earthquakes > 0))])
        FP = len(stress[((stress>=CFSth) & (earthquakes == 0))])
        FN = len(stress[((stress <CFSth) & (earthquakes > 0))])
        TN = len(stress[((stress <CFSth) & (earthquakes == 0))])
#        print('TP :',TP)
#        print('FP :' ,FP)
#        print('FN :',FN)
#        print('TN :',TN)
        if FN > 0 or TP > 0:
            TPR = TP/float(TP + FN)
            FPR = FP/float(TN + FP)
            ROCdata.append([FPR, TPR, CFSth])  #[x-axis data, y-axis data, threshold]
        if FN > 0 and TN > 0:
            MCC = (((TP*TN)-(FP*FN))/((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))**0.5 +1)/2 #MCC rescaled between 0-1
#            print(MCC)
            F1 = (2*TP)/((2*TP)+FP+FN)       #F1-measure
            MCC_F1.append([F1, MCC, CFSth])    #[x-axis data, y-axis data, threshold]    
    
    return numpy.array(ROCdata), numpy.array(MCC_F1)


def calc_mean_dist(x1, y1, x2, y2):
    return numpy.mean(numpy.sqrt((x1-x2)**2 + (y1-y2)**2))

def calc_mccf1_metric(f1,mcc, min_dist=False):
    """
    Calculate mcc-f1  metric
    f1: list or array of f1 metric (y-axis)
    mcc: list or array of mcc (x-axis)
    """
    if min_dist:
        #Use min distance to compute metric. 
        dist = numpy.sqrt( (f1-1)**2 + (mcc-1)**2 )
        MCC_F1_metric = 1 - min(dist)/numpy.sqrt(2)  #2 
    else:
        #Use mean-based metric defined in the paper
        max_mc_index = numpy.argmax(mcc)
            
        # divide the MCC-F1 curve in two parts
        # before the max MCC and after the max MCC
        
        # ---------- before
        mcc_before = mcc[:max_mc_index]
        f_before = f1[:max_mc_index]
        # ---------- after
        mcc_after = mcc[max_mc_index:]
        f_after = f1[max_mc_index:]
            
        #get length of bin
        bins = 1000
        unit_len = (max(mcc) - min(mcc))/bins
        
        #Now in each bin (remember these bins are vertical) calculate the mean distance of MCC-F1
        # points to the point of perfection (1,1). Then, take average of all mean points
        mean_dist_before, mean_dist_after = [], []
        for i in range (bins) :
        #Bin-wise average of before
            mcc_indx_before = ((mcc_before >= min(mcc) + i*unit_len) & (mcc_before <= min(mcc) + (i+1)*unit_len))
            mean_dist_before = numpy.append(mean_dist_before, calc_mean_dist(f_before[mcc_indx_before], mcc_before[mcc_indx_before], 1, 1))
                
            #Bin-wise average after 
            mcc_indx_after = ((mcc_after >= min(mcc) + i*unit_len) & (mcc_after <= min(mcc) + (i+1)*unit_len))
            mean_dist_after = numpy.append(mean_dist_after, calc_mean_dist(f_after[mcc_indx_after], mcc_after[mcc_indx_after], 1, 1))
            
        # club all distance in one array
        all_mean_distances = numpy.concatenate((mean_dist_before,mean_dist_after))
    
        # drop NaN
        all_mean_distances_no_NaN = all_mean_distances[~numpy.isnan(all_mean_distances)]
        
        # take average
        distance_average = numpy.mean(all_mean_distances_no_NaN)
        
        # calculate MCC-F1 metric
        MCC_F1_metric = 1 - distance_average/numpy.sqrt(2)  #2
        
    return MCC_F1_metric
