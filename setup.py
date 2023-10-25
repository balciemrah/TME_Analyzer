import cx_Freeze
import sys
import os
# import tkinter
# import tkinter.messagebox
# import tkinter.filedialog
# import skimage.io
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib
# import matplotlib.backends.backend_tkagg
# from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as NT2Tk
import scipy.ndimage
import xlwt
# from math import atan2, pi as PI
import scipy
import scipy.spatial
# import skimage
# import skimage.io
# import skimage.measure
# import skimage.feature
# import skimage.morphology.watershed
# import skimage.measure._moments
# import openTSNE
# import multiprocessing
# import multiprocessing.pool
# import numba
# import llvmlite
from os import path
# import google.protobuf
# import tensorboard

# import stardist
# import csbdeep
# import tensorflow

base = None

if sys.platform == 'win32':
    base = "Win32GUI"

executables = [cx_Freeze.Executable("TME_analyzer_test.py", base=base)]

#os.environ['TCL_LIBRARY'] = r'C:\Users\TME facility\AppData\Local\Programs\Python\Python38\tcl\tcl8.6'
#os.environ['TK_LIBRARY'] = r'C:\Users\TME facility\AppData\Local\Programs\Python\Python38\tcl\tk8.6'


includefiles_list=[]
# skimage_path = path.dirname(skimage.__file__)
# includefiles_list.append(path.dirname(scipy.__file__))
# includefiles_list.append(skimage_path)
# multiprocessing_path = path.dirname(multiprocessing.__file__)
# includefiles_list.append(multiprocessing_path)
# includefiles_list.append(path.dirname(xlwt.__file__))
# numba_path = path.dirname(numba.__file__)
# includefiles_list.append(numba_path)
# llvmlite_path = path.dirname(llvmlite.__file__)
# includefiles_list.append(llvmlite_path)
# openTSNE_path = path.dirname(openTSNE.__file__)
# includefiles_list.append(openTSNE_path)
# tensorflow_path = path.dirname(tensorflow.__file__)
# includefiles_list.append(tensorflow_path)
# tensorboard_path = path.dirname(tensorboard.summary._tf.summary.__file__)
# includefiles_list.append(tensorboard_path)

cx_Freeze.setup(
    name = "Image Analysis",
    options = {"build_exe": {"packages":["tensorboard.summary._tf.summary","tensorboard","google","tensorflow","stardist","csbdeep","llvmlite","numba", "tkinter", "numpy", "matplotlib", "skimage", "scipy", "pandas", "skimage.io", "skimage.measure", "skimage.feature", "skimage.morphology", "matplotlib.pyplot","os"], "includes":["matplotlib.backends.backend_tkagg"], "include_files": includefiles_list}},
    version = "0.01",
    description = "first trial",
    executables = executables
    )