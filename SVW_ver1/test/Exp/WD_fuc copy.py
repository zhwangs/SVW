# here

'''
You need to install the following packages:

pip install numpy
pip install scipy
pip install h5py
pip install pandas

pip install matplotlib


pip install pypandoc

'''
# %% importing modules
#################################################################
#################################################################

# general ------------------------------------------------------------
import numpy as np
import scipy as sp
import h5py
import re
import time
import datetime
import os
import pathlib
from numpy.polynomial.legendre import legval, legder
from scipy.special import lpmv, lpmn, clpmn, jv, jvp, yv, yvp, hankel1, hankel2, h1vp, h2vp, factorial, sph_harm
import sys
import pandas as pd

from pathlib import Path

import json

import runpy

from time import strftime,time
from scipy.special import spherical_jn, spherical_yn

# plotting -------------------------------------------------------------
import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

font = {    
        'size'   : 5}

matplotlib.rc('font', **font)

import torch 
from scipy.special import factorial
from torch.autograd import grad


#%% Rotation 

def rotor(x,y,z,axis)