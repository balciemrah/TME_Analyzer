import cx_Freeze
import sys
import os
import scipy.ndimage
import xlwt
import scipy
import scipy.spatial
from os import path

base = None

if sys.platform == 'win32':
    base = "Win32GUI"

executables = [cx_Freeze.Executable("TME_analyzer_test.py", base=base)]

cx_Freeze.setup(
    name = "Image Analysis",
    options = {"build_exe": {"packages":["tensorboard.summary._tf.summary","tensorboard","google",
                                         "tensorflow","stardist","csbdeep","llvmlite","numba", "tkinter", 
                                         "numpy", "matplotlib", "skimage", "scipy", "pandas", "skimage.io",
                                        "skimage.measure", "skimage.feature", "skimage.morphology",
                                          "matplotlib.pyplot","os"], "includes":["matplotlib.backends.backend_tkagg"]}},
    version = "0",
    description = "test-pack",
    executables = executables
    )