# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 14:41:55 2021

@author: Furkan
"""


import importlib.machinery

pyMCDS = importlib.machinery.SourceFileLoader('pyMCDS','./analysis/pyMCDS.py').load_module()

import os
from pathlib import Path
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
from fury import window, actor, utils, primitive, io, ui
from fury.data import read_viz_textures, fetch_viz_textures
import itertools
import vtk
import glob
import time
import random
import scipy.io as sio




def plot_micenvs (time_point):
    fine_data = sio.loadmat(time_point + "_microenvironment0.mat")['multiscale_microenvironment']
    # coarse_data = sio.loadmat(time_point + "_microenvironment1.mat")['multiscale_microenvironment']
    dx = fine_data[0,1]-fine_data[0,0]
    fine_x = np.unique(fine_data[0,:])
    fine_y = np.unique(fine_data[1,:])
    fine_X, fine_Y = np.meshgrid(fine_x, fine_y)
    fine_oxy = fine_data[4,np.where(fine_data[2,:] == 16)]
    fine_oxy = fine_oxy.reshape((len(fine_y),len(fine_x)))

    fig,axs = plt.subplots(1,1)
    cp = axs.contourf(fine_X,fine_Y,fine_oxy)
    axs.axis('equal')
    axs.set(xlim=(-3000,3000), ylim=(-500,500))
    fig.colorbar(cp,format = "%f")
    fig.tight_layout()
    plt.pause(0.5)
    return fine_data





main_path = Path(os.getcwd()).parent
out_path = os.path.join(main_path, "output")

os.chdir(out_path)

time_point = "output00000002"
data = plot_micenvs(time_point)


data1 = data[4,np.where(data[0,:]) == -2864]
print(data1)

# for i in range(0,10):
#     timepoint=time_point+str(i)
#     plot_micenvs(timepoint)
    
