# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 20:05:18 2020

@author: Patrick
"""
import numpy as np
from matplotlib import pyplot as plt

frequencies = np.array([0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.15])
omegas = frequencies * 2* np.pi
amplitudes = np.array([])