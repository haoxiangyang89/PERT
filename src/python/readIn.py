#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 15:20:53 2019

@author: haoxiangyang
"""
import os
import pickle
import csv
import numpy as np
from defClass import *

distrList = ['normal','exponential','lognormal','uniform']

# read in the project activity information and link information
def readInP(pInputAdd,kInputAdd):
    fi = open(pInputAdd,"r")
    csvReader = csv.reader(fi,dialect = 'excel')
    pRaw = []
    for item in csvReader:
        pRaw.append(item)
    fi.close()
    
    fi = open(kInputAdd,"r")
    csvReader = csv.reader(fi,dialect = 'excel')
    kRaw = []
    for item in csvReader:
        kRaw.append(item)
    fi.close()
    
    np = len(pRaw)
    mp = len(pRaw[0])
    nk = len(kRaw)
    
    # total budget and nominal scenario probability
    B = float(pRaw[0][0])
    p0 = float(pRaw[0][1])
    
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
    for i in range(1,np):
        lineID = int(pRaw[i][0])
        II.append(lineID)
        D[lineID] = float(pRaw[i][1])
        
        jStart = 2
        Ji[lineID] = []
        b[lineID] = {}
        eff[lineID] = {}
        jCounter = 0
        while jStart <= mp-1:
            jCounter += 1
            Ji[lineID].append(jCounter)
            b[lineID][jCounter] = float(pRaw[i][2*jCounter])
            eff[lineID][jCounter] = float(pRaw[i][1+2*jCounter])
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
        K.append((fromI,toI))
        Pre[toI].append(fromI)
        Succ[fromI].append(toI)
        
    pData = pInfo(II,Ji,D,b,eff,B,p0,K,Pre,Succ)
    return pData

def readInUnc(phiInputAdd):
    fi = open(phiInputAdd,"r")
    csvReader = csv.reader(fi,dialect = 'excel')
    phiRaw = []
    for item in csvReader:
        phiRaw.append(item)
    fi.close()
    
    phin = len(phiRaw)
    nameD = phiRaw[0][0]
    dparams = {}
    for n in range(1,phin):
        dparams[int(phiRaw[n][0])] = [float(i) for i in phiRaw[n][1:] if i != '']
    return nameD,dparams

def readInH(hInputAdd):
    fi = open(hInputAdd,"r")
    csvReader = csv.reader(fi,dialect = 'excel')
    hRaw = []
    for item in csvReader:
        hRaw.append(item)
    fi.close()
    nameH = hRaw[0][0]
    Hparams = []
    for i in range(len(hRaw[1])):
        if hRaw[1][i] != '':
            Hparams.append(float(hRaw[1][i]))
    return nameH,Hparams

def ranGen(nameR,Rparams,size = 0):
    nameR = nameR.lower()
    if nameR in distrList:
        evalString = "np.random.{}(".format(nameR)
        for ritem in Rparams:
            evalString += "{},".format(ritem)
        if size == 0:
            evalString = evalString[:-1]+')'
            returnedVal = eval(evalString)
        else:
            evalString = evalString[:-1]+'{})'.format(size)
            returnedVal = eval(evalString)
    else:
        if nameR == 'singleton':
            returnedVal = Rparams
    return returnedVal

def autoUGen(nameH, Hparams, nameD, dparams, Omegan, totalProb):
    disData = {}
    Omega = list(range(Omegan))
    for omega in Omega:
        H = round(ranGen(nameH,Hparams),4)
        d = {}
        for i in dparams.keys():
            d[i] = round(ranGen(nameD,dparams[i]),4)
        pomega = totalProb/Omegan
        disData[omega] = disInfo(H,d,pomega)
    return disData,Omega

def orderdisData(disData,Omega):
    dHList = [disData[omega].H for omega in Omega]
    omegaOrdered = sorted(range(len(dHList)), key=lambda k: dHList[k])
    disDataNew = {}
    for omega in Omega:
        disDataNew[omega] = disData[omegaOrdered[omega]]
    return disDataNew

def genData(filePath,Omegasize,dataSize = 1,pName = 'test_P.csv',kName = 'test_K.csv',\
            phiName = 'test_Phi.csv',hName = 'test_H.csv',dOnly = 0,hOnly = 0,saveOpt = 0):
    pInputAdd = os.path.join(filePath,pName)
    kInputAdd = os.path.join(filePath,kName)
    phiInputAdd = os.path.join(filePath,phiName)
    hInputAdd = os.path.join(filePath,hName)
    
    pData = readInP(pInputAdd,kInputAdd)
    nameD,dparams = readInUnc(phiInputAdd)
    nameH,Hparams = readInH(hInputAdd)
    
    disDataSet = []
    for ds in range(dataSize):
        if (dOnly == 0)and(hOnly == 0):
            disData,Omega = autoUGen(nameH,Hparams,nameD,dparams,Omegasize,1 - pData.p0)
            disData = orderdisData(disData,Omega)
        elif (dOnly != 0)and(hOnly == 0):
            distrD = {}
            for i in pData.II:
                if nameD == "LogNormal":
                    distrD[i] = np.exp(dparams[0]+dparams[1]**2/2)
                elif (nameD == "Normal")or(nameD == "Exponential"):
                    distrD[i] = dparams[0]
                elif nameD == "Uniform":
                    distrD[i] = (dparams[0] + dparams[1])/2
            disData,Omega = autoUGen(nameH,Hparams,"singleton",distrD,Omegasize,1 - pData.p0)
        elif (dOnly == 0)and(hOnly != 0):
            if nameH == "LogNormal":
                distrH = np.exp(Hparams[0]+ Hparams[1]**2/2)
            elif (nameH == "Normal")or(nameH == "Exponential"):
                distrH = Hparams[0]
            elif nameH == "Uniform":
                distrH = (Hparams[0] + Hparams[1])/2
            disData,Omega = autoUGen("singleton",distrH,nameD,dparams,Omegasize,1 - pData.p0)
        disDataSet.append(disData)
    if saveOpt == 0:
        return pData,disDataSet,nameD,nameH,dparams,Hparams
    else:
        with open(os.path.join(filePath,'probScen.p'), 'wb') as fp:
            pickle.dump([pData,disDataSet], fp, protocol=pickle.HIGHEST_PROTOCOL)