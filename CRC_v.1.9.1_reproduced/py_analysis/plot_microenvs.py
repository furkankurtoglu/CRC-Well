# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 14:41:55 2021

@author: Furkan
"""


import importlib.machinery

pyMCDS = importlib.machinery.SourceFileLoader('pyMCDS','./analysis/pyMCDS.py').load_module()


import numpy as np
import pandas as pd
from fury import window, actor, utils, primitive, io, ui
from fury.data import read_viz_textures, fetch_viz_textures
import itertools
import vtk
import glob
import time
import random




output_folder = "..\\output\\"
outpur_file = "" #"output00000000.xml"
output = output_folder + outpur_file

mcds=pyMCDS.pyMCDS(output)
data=mcds.get_cell_df()