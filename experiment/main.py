#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 18:54:26 2023

@author: khawajasim
"""
import os
import experiment.square_circle_grid_tests as single_run
import experiment.square_circle_grid_tests_multi_run as multi_run
import experiment.aftershocks_tests_real_data as real_test
import experiment.analyze_min_required_data as min_data

folder = '../output/'

if not os.path.exists(folder):
    os.mkdir(folder)

def run():
    single_run.main()
    multi_run.main()
    min_data.main()
    real_test.main()

