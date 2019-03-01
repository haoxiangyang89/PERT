#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 17:45:06 2019

@author: haoxiangyang
"""

# functions to tighten the partition

from defClass import *
from gurobipy import *
from tighten import *
import copy
import numpy as np
from extForm import *
from ubCalFunc import *

def findSuccAll(pData):
    # find all successors of activity i (children, children's children)
    allSucc = {}
    for i in pData.II:
        succList = pData.Succ[i].copy()
        succAll = []
        while succList != []:
            currentSucc = succList[0]
            succAll.append(currentSucc)
            succList.pop(0)
            for j in pData.Succ[currentSucc]:
                if not((j in succAll) or (j in succList)):
                    succList.append(j)
        allSucc[i] = succAll
    return allSucc

def detCal(pData,ii,jj):
    # find the distance between i and j under nominal cases without crashing
    mp = Model('detCal')
    mp.Params.OutputFlag = 0
    t = mp.addVars(pData.II, lb = 0)
    JJ = [(i,j) for i in pData.II for j in pData.Ji[i]]
    x = mp.addVars(JJ, lb = 0, ub = 1)
    
    budgetConstr = mp.addConstr(sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B)
    durationConstr = mp.addConstrs(t[k[1]] - t[k[0]] >= pData.D[k[0]]*(1-sum(pData.eff[k[0]][j]*x[k[0],j] for j in pData.Ji[k[0]])) for k in pData.K)
    xConstr = mp.addConstrs(sum(x[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II)
    
    mp.setObjective(t[jj] - t[ii], GRB.MINIMIZE)
    mp.update()
    mp.optimize()
    disIJ = mp.objVal
    
    return disIJ
   
def revPar(prevPartSet, currentPart):
    # i specific
    corrPart = 0
    for par in range(len(prevPartSet)):
        if (currentPart.startH >= prevPartSet[par].startH)and(currentPart.endH <= prevPartSet[par].endH):
            corrPart = par
    return corrPart

def iniPart(pData,disData,Omega,sN,MM,returnOpt = 0):
    Tmax = disData[len(Omega) - 1].H + longestPath(pData)[0]
    H = {}
    H[0] = 0
    H[len(Omega) + 1] = Tmax
    for omega in Omega:
        H[omega + 1] = disData[omega].H
    textList = []
    xextList = []

    # sample sN scenarios and solve the small extensive formulation
    tList = {}
    for i in pData.II:
        tList[i] = []
    ubList = []
    ubMin = np.inf
    tBest = {}
    xBest = {}
    θBest = {}
    for m in range(MM):
        #randList = rand(Omega,sN)
        disData1 = {}
        Omega1 = range(sN)
        randList = [omega*MM + m for omega in Omega1]
        for i in Omega1:
            disData1[i] = copy.deepcopy(disData[randList[i]])
            disData1[i].prDis = (1 - pData.p0)/len(Omega1)
        text,xext,fext,gext,mext = extForm(pData,disData1,Omega1,1e-4,999999)
        ubext,comegaList = ubCal(pData,disData,Omega,xext,text,999999,1)
        ubList.append(ubext)
        textList.append(text)
        xextList.append(xext)
        # record the best solution
        if ubext < ubMin:
            ubMin = ubext
            for i in pData.II:
                tBest[i] = text[i]
                for j in pData.Ji[i]:
                    xBest[i,j] = xext[i,j]
            for omega in Omega:
                θBest[omega] = comegaList[omega]
        for i in pData.II:
            tList[i].append(text[i])

    # return the min-max decomposition
    tHList = []
    for i in pData.II:
        tHstart = max([omega for omega in H.keys() if H[omega] <= min(tList[i])])
        tHend = min([omega for omega in H.keys() if H[omega] >= max(tList[i])])
        tHList.append([i,tHstart,tHend])
    
    if returnOpt == 0:    
        return ubList,tHList,ubMin,tBest,xBest,θBest
    else:
        return ubList,tHList,ubMin,tBest,xBest,θBest,textList,xextList

def splitPar_CI(PartSet,PartDet,splitInfo):
    # i generic
    newPartSet = copy.deepcopy(PartSet)
    newPartDet = copy.deepcopy(PartDet)
    for (i,splitStart,splitEnd) in splitInfo:
        partSetiTemp = []
        partDetiTemp = []
        for par in range(len(PartSet[i])):
            # split the current partition, essentially split twice using the start and the end
            if (splitStart < PartSet[i][par].endH)and(splitStart > PartSet[i][par].startH)and(PartDet[i][par] == 0):
                if (splitEnd < PartSet[i][par].endH)and(splitEnd > PartSet[i][par].startH):
                    part1 = partType(PartSet[i][par].startH,splitStart)
                    part2 = partType(splitStart,splitEnd)
                    part3 = partType(splitEnd,PartSet[i][par].endH)
                    partSetiTemp.append(part1)
                    partSetiTemp.append(part2)
                    partSetiTemp.append(part3)
                    partDetiTemp.append(0)
                    partDetiTemp.append(0)
                    partDetiTemp.append(0)
                else:
                    part1 = partType(PartSet[i][par].startH,splitStart)
                    part2 = partType(splitStart,PartSet[i][par].endH)
                    partSetiTemp.append(part1)
                    partSetiTemp.append(part2)
                    partDetiTemp.append(0)
                    partDetiTemp.append(0)
            else:
                if (splitEnd < PartSet[i][par].endH)and(splitEnd > PartSet[i][par].startH):
                    part1 = partType(PartSet[i][par].startH,splitEnd)
                    part2 = partType(splitEnd,PartSet[i][par].endH)
                    partSetiTemp.append(part1)
                    partSetiTemp.append(part2)
                    partDetiTemp.append(0)
                    partDetiTemp.append(0)
                else:
                    partSetiTemp.append(PartSet[i][par])
                    partDetiTemp.append(PartDet[i][par])
        newPartSet[i] = partSetiTemp
        newPartDet[i] = partDetiTemp
    return newPartSet,newPartDet

def splitPar3(PartSet,PartDet,splitInfo):
    # split the current fractional partition into equal 3 parts
    # i generic
    newPartSet = copy.deepcopy(PartSet)
    newPartDet = copy.deepcopy(PartDet)
    for (i,splitStart,splitEnd) in splitInfo:
        partSetiTemp = []
        partDetiTemp = []
        for par in range(len(PartSet[i])):
            # split the current partition
            if (splitEnd <= PartSet[i][par].endH)and(splitStart >= PartSet[i][par].startH):
                mod3 = (PartSet[i][par].endH - PartSet[i][par].startH) % 3
                len3 = (PartSet[i][par].endH - PartSet[i][par].startH) // 3
                if len3 >= 1:
                    # if the partition is larger than 3 items in between
                    if mod3 == 0:
                        part1 = partType(PartSet[i][par].startH,PartSet[i][par].startH + len3)
                        part2 = partType(PartSet[i][par].startH + len3,PartSet[i][par].startH + len3*2)
                        part3 = partType(PartSet[i][par].startH + len3*2,PartSet[i][par].endH)
                    elif mod3 == 1:
                        part1 = partType(PartSet[i][par].startH,PartSet[i][par].startH + len3)
                        part2 = partType(PartSet[i][par].startH + len3,PartSet[i][par].startH + len3*2)
                        part3 = partType(PartSet[i][par].startH + len3*2,PartSet[i][par].endH)
                    else:
                        part1 = partType(PartSet[i][par].startH,PartSet[i][par].startH + len3)
                        part2 = partType(PartSet[i][par].startH + len3,PartSet[i][par].startH + len3*2 + 1)
                        part3 = partType(PartSet[i][par].startH + len3*2 + 1,PartSet[i][par].endH)
                    partSetiTemp.append(part1)
                    partSetiTemp.append(part2)
                    partSetiTemp.append(part3)
                    partDetiTemp.append(0)
                    partDetiTemp.append(0)
                    partDetiTemp.append(0)
                else:
                    splitPos = int(np.floor((splitStart + splitEnd)/2))
                    part1 = partType(PartSet[i][par].startH,splitPos)
                    part2 = partType(splitPos,PartSet[i][par].endH)
                    partSetiTemp.append(part1)
                    partSetiTemp.append(part2)
                    partDetiTemp.append(0)
                    partDetiTemp.append(0)
            else:
                partSetiTemp.append(PartSet[i][par])
                partDetiTemp.append(PartDet[i][par])
        newPartSet[i] = partSetiTemp
        newPartDet[i] = partDetiTemp
    return newPartSet,newPartDet