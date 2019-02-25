#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  2 17:16:49 2019

@author: haoxiangyang
"""

# test the python script for PERT
import csv
import numpy as np
import os
import pickle
import copy

os.chdir('/Users/haoxiangyang/Desktop/Git/PERT/src/python/')
from defClass import *
from readIn import *
from extForm import *
from ubCalFunc import *
from tighten import *

OmegaSize = 1
pData,disDataSet,nameD,nameH,dparams,Hparams = genData('/Users/haoxiangyang/Desktop/PERT_tests/14_Lognormal_Exponential/',OmegaSize)
Omega = list(range(OmegaSize))
disData = disDataSet[0]

#%%

pdData = copy.deepcopy(pData)
for i in pData.II:
    if i != 0:
        pdData.D[i] = pData.D[i] + max([disData[omega].d[i] for omega in Omega])
    else:
        pdData.D[i] = pData.D[i]
        
lDict = longestPath(pdData)
Tmax1 = lDict[0]

tdet,xdet,fdet = detBuild(pData)
ubdet = ubCal(pData,disData,Omega,xdet,tdet,Tmax1)

eH = np.mean([disData[omega].H for omega in Omega])
ed = {}
for i in pData.II:
    if i != 0:
        ed[i] = np.mean([disData[omega].d[i] for omega in Omega])
texp,xexp,fexp,Gexp,mexp = expModel(pData,eH,ed,Tmax1)
ubexp = ubCal(pData,disData,Omega,xexp,texp,Tmax1)

disData1 = copy.deepcopy(disData)
for omega in Omega:
    disData1[omega].H = eH
tdOnly,xdOnly,fdOnly,gdOnly,mdOnly = extForm(pData,disData1,Omega,1e-4,Tmax1)
ubdOnly = ubCal(pData,disData,Omega,xdOnly,tdOnly,Tmax1)

disData1 = copy.deepcopy(disData)
for omega in Omega:
    disData1[omega].d = ed
tHOnly,xHOnly,fHOnly,gHOnly,mHOnly = extForm(pData,disData1,Omega,1e-4,Tmax1)
ubHOnly = ubCal(pData,disData,Omega,xHOnly,tHOnly,Tmax1)

tFull,xFull,fFull,gFull,mFull = extForm(pData,disData,Omega,1e-4,Tmax1);
ubFull = ubCal(pData,disData,Omega,xFull,tFull,Tmax1);