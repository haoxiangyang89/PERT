#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 22:54:08 2018

@author: haoxiangyang
"""

# read in the data
from defType import *
import numpy as np
import scipy.stats
import csv
import pickle
import os

def buildDistrn(nameDistr,paramDistr):
    if nameDistr == "Exponential":
        μ = paramDistr[0]
        distrObj = scipy.stats.expon(0,μ)
    elif nameDistr == "LogNormal":
        μ = paramDistr[0]
        σ = paramDistr[1]
        distrObj = scipy.stats.lognorm([σ],scale = np.exp(μ))
    elif nameDistr == "Gamma":
        α = paramDistr[0]
        θ = paramDistr[1]
        distrObj = scipy.stats.gamma(scale = α,loc = θ)
    elif nameDistr == "Uniform":
        la = paramDistr[0]
        ub = paramDistr[1]
        distrObj = scipy.stats.uniform(loc = la,scale = ub+1e-10 - la)
    elif nameDistr == "Normal":
        μ = paramDistr[0]
        σ = paramDistr[1]
        distrObj = scipy.stats.norm(loc = μ,scale = σ)
    elif nameDistr == "Singleton":
        la = paramDistr[0]
        distrObj = scipy.stats.uniform(loc = la,scale = 0)
        
    return distrObj

# input function that reads in the project data
def readInP(pInputAdd,kInputAdd):
    fi = open(pInputAdd,"r")
    csvReader = csv.reader(fi)
    pRaw = []
    for item in csvReader:
        pRaw.append(item)
    fi.close()
    
    fi = open(kInputAdd,"r")
    csvReader = csv.reader(fi)
    kRaw = []
    for item in csvReader:
        kRaw.append(item)
    fi.close()
    
    npp = len(pRaw)
    nk = len(kRaw)
    
    # total budget
    B = float(pRaw[0][0])
    # nominal scenario probability
    p0 = float(pRaw[0][1])
    
    # 0 as the last activity
    II = [0]
    Ji = {}
    D = {}
    D[0] = 0
    Ji[0] = []
    
    b = {}
    b[0] = {}
    eff = {}
    eff[0] = {}
    
    # read in the activity information
    for i in range(1,npp):
        lineID = int(pRaw[i][0])
        II.append(lineID)
        D[lineID] = float(pRaw[i][1])
        
        jStart = 2
        Ji[lineID] = []
        b[lineID] = {}
        eff[lineID] = {}
        jCounter = 1
        while jStart <= len(pRaw[i])-1:
            Ji[lineID].append(jCounter)
            b[lineID][jCounter] = float(pRaw[i][jStart])
            eff[lineID][jCounter] = float(pRaw[i][1+jStart])
            jCounter += 1
            jStart += 2
            
    # read in the activity precedence information
    K = []
    Pre = {}
    Pre[0] = []
    Succ = {}
    Succ[0] = []
    for i in II:
        if i != 0:
            K.append((i,0))
            Pre[0].append(i)
            Pre[i] = []
            Succ[i] = [0]
            
    for k in range(nk):
        fromI = int(kRaw[k][0])
        toI = int(kRaw[k][1])
        K.append((fromI, toI))
        Pre[toI].append(fromI)
        Succ[fromI].append(toI)
    
    pData = pInfo(II,Ji,D,b,eff,B,p0,K,Pre,Succ)
    return pData

# read in the scenario information from file
def readInDis(OmegaInputAdd):
    fi = open(OmegaInputAdd,"r")
    csvReader = csv.reader(fi)
    OmegaRaw = []
    for item in csvReader:
        OmegaRaw.append(item)
    fi.close()
    
    noo = len(OmegaRaw)
    moo = len(OmegaRaw[0])
    Omega = [i for i in range(1,moo)]
    disData = {}
    for omega in Omega:
        d = {}
        H = float(OmegaRaw[0][omega])
        pomega = float(OmegaRaw[1][omega])
        for n in range(2,noo):
            d[int(OmegaRaw[n][0])] = float(OmegaRaw[0][omega])
            
        disData[omega] = disInfo(H,d,pomega)
        
    return disData,Omega

# read in the magnitude distribution information from file
def readInUnc(phiInputAdd):
    fi = open(phiInputAdd,"r")
    csvReader = csv.reader(fi)
    phiRaw = []
    for item in csvReader:
        phiRaw.append(item)
    fi.close()
    
    phin = len(phiRaw)
    phim = len(phiRaw[0])
    
    nameD = phiRaw[0][0]
    dparams = {}
    for n in range(1,phin):
        dparams[int(phiRaw[n][0])] = [float(phiRaw[n][i]) for i in range(1,phim) if phiRaw[n][i] != '']
    return nameD,dparams

# read in the timing distribution information from file
def readH(hInputAdd):
    fi = open(hInputAdd,"r")
    csvReader = csv.reader(fi)
    hRaw = []
    for item in csvReader:
        hRaw.append(item)
    fi.close()
    nameH = hRaw[0][0]
    mh = len(hRaw[0])
    Hparams = []
    for i in range(mh):
        Hparams.append(float(hRaw[1][i]))
    return nameH,Hparams

# automatically uncertainty data generation
def autoUGen(nameH, Hparams, nameD, dparams, Omegan, totalProb):
    disData = {}
    Omega = [i for i in range(Omegan)]
    distrnH = buildDistrn(nameH,Hparams);
    
    distrnD = {}
    for i in dparams.keys():
        distrnD[i] = buildDistrn(nameD,dparams[i])
        
    for omega in Omega:
        H = round(distrnH.rvs(),4)
        d = {}
        for i in dparams.keys():
            d[i] = round(distrnD[i].rvs(),4)
        pomega = totalProb/Omegan
        disData[omega] = disInfo(H,d,pomega)
        
    return disData,Omega

# order the disData by the timing
def orderdisData(disData,Omega):
    dHList = [disData[omega].H for omega in Omega]
    omegaOrdered = np.argsort(dHList)
    disDataNew = {}
    for omega in range(len(Omega)):
        disDataNew[omega] = disData[omegaOrdered[omega]]
        
    return disDataNew

# generate the data in a folder
def genData(filePath,Omegasize,dataSize = 1,pName = "test_P.csv",kName = "test_K.csv",\
            phiName = "test_Phi.csv",hName = "test_H.csv", dOnly = 0, hOnly = 0, saveOpt = 0):
    pInputAdd = os.path.join(filePath,pName)
    kInputAdd = os.path.join(filePath,kName)
    phiInputAdd = os.path.join(filePath,phiName)
    hInputAdd = os.path.join(filePath,hName)
    
    pData = readInP(pInputAdd,kInputAdd)
    nameD,dparams = readInUnc(phiInputAdd)
    nameH,Hparams = readH(hInputAdd)
    
    disDataSet = []
    
    for ds in range(dataSize):
        if (dOnly == 0)and(hOnly == 0):
            disData,Omega = autoUGen(nameH,Hparams,nameD,dparams,Omegasize,1 - pData.p0)
            disData = orderdisData(disData,Omega)
        elif (dOnly != 0)and(hOnly == 0):
            distrD = {}
            for i in pData.II:
                distrnDtemp = buildDistrn(nameD,dparams[i])
                distrD[i] = distrnDtemp.mean()
            disData,Omega = autoUGen(nameH,Hparams,"Singleton",distrD,Omegasize,1 - pData.p0)
            disData = orderdisData(disData,Omega)
        elif (dOnly == 0)and(hOnly != 0):
            distrnHtemp = buildDistrn(nameH,Hparams)
            distrH = distrnHtemp.mean()
            disData,Omega = autoUGen("Singleton",distrH,nameD,dparams,Omegasize,1 - pData.p0)
            disData = orderdisData(disData,Omega)
        disDataSet.append(disData)
    if saveOpt == 0:
        return pData,disDataSet,nameD,nameH,dparams,Hparams
    else:
        dataDump = {"pData":pData,"disDataSet":disDataSet}
        picklePath = os.path.join(filePath,"test.p")
        with open(picklePath, 'wb') as fp:
            pickle.dump(dataDump, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
# load the pickle data
def loadData(filePath, fName = "test.p"):
    fileAdd = os.path.join(filePath, fName)
    disDataSetR = pickle.load(open(fileAdd, 'rb'))
    disDataSet = disDataSetR["disDataSet"]
    return disDataSet