#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 18:54:26 2023

@author: khawajasim
"""
import os
import square_circle_grid_tests, square_circle_grid_tests_multi_run, aftershocks_tests_real_data

folder = '../output/'

if not os.path.exists(folder):
    os.mkdir(folder)

def run():
    square_circle_grid_tests.main()
    square_circle_grid_tests_multi_run.main()
    aftershocks_tests_real_data.main()
