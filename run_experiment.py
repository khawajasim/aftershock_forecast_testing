#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 22:15:54 2023

@author: khawajasim
"""

# Local modules
# import sys
from experiment import run
import time


if __name__ == "__main__":
    start_time = time.time()
    run()
    print("Time to Generate all Figures (mins) :" % (time.time() - start_time)/60)


