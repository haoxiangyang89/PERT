#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:10:33 2019

@author: haoxiangyang
"""

# auxiliary functions to tighten the branch-and-cut process

import numpy as np
from gurobipy import *

def longestPath(pData):
    # obtain the longest path using negative Dijkstra
    lDict = {}
    finishedList = []
    activeList = []
    for i in pData.II:
        if pData.Pre[i] == []:
            lDict[i] = 0
        else:
            lDict[i] = -np.inf
    while len(finishedList) < len(pData.II):
        for i in pData.II:
            if not((i in activeList) or (i in finishedList)):
                enterBool = True
                for j in pData.Pre[i]:
                    if not(j in finishedList):
                        enterBool = False
                if enterBool:
                    activeList.append(i)
        # find the largest activity in activeList
        maxVal = -np.inf
        maxInd = -1
        for iInd in activeList:
            if maxVal < lDict[iInd]:
                maxVal = lDict[iInd]
                maxInd = iInd
        # update the activities connected to maxInd
        for j in pData.Succ[maxInd]:
            if lDict[j] < lDict[maxInd] + pData.D[maxInd]:
                lDict[j] = lDict[maxInd] + pData.D[maxInd]
        finishedList.append(maxInd)
        activeList.remove(maxInd)
    return lDict

def iSolve(pData,iTarget,elInd,ubTemp = 9999999):
    # build the LP to obtain the earliest/latest possible starting time of an activity
    mp = Model('iSolve')
    mp.Params.OutputFlag = 0
    
    JJ = [(i,j) for i in pData.II for j in pData.Ji[i]]
    t = mp.addVars(pData.II, lb = 0, name = "t")
    x = mp.addVars(JJ, lb = 0, ub = 1, name = "x")
    
    durationConstr = mp.addConstrs((t[k[1]] - t[k[0]] >= pData.D[k[0]]*(1 - sum(pData.eff[k[0]][j]*x[k[0],j] for j in pData.Ji[k[0]])) for k in pData.K), name = "durationConstr")
    budgetConstr = mp.addConstr((sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B), name = "budgetConstr")
    xConstr = mp.addConstrs((sum(x[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II), name = "xConstr")
    ubCon = mp.addConstr((t[0] <= ubTemp), name = "ubCon")
    
    if elInd == 0:
        mp.setObjective(t[iTarget], GRB.MINIMIZE)
    else:
        mp.setObjective(t[iTarget], GRB.MAXIMIZE)

    mp.update()
    mp.optimize()
    return mp.objVal

def precludeRel(pData,H,ub = np.inf):
    # solve |II| number of LP to obtain each activity's earliest starting time
    # initialize the brInfo
    brInfo = np.zeros([len(pData.II), len(H) - 2])
    for i in pData.II:
        tearly = iSolve(pData,i,0,ub)
        for hIter in H.keys():
            if (hIter != 0)and(hIter != max(H.keys())):
                if tearly >= H[hIter]:
                    brInfo[pData.II.index(i),hIter - 1] = 1
    if ub != np.inf:
        for i in pData.II:
            tlate = iSolve(pData,i,1,ub)
            for hIter in H.keys():
                if (hIter != 0)and(hIter != max(H.keys())):
                    if tlate <= H[hIter]:
                        brInfo[pData.II.index(i),hIter - 1] = -1
    return brInfo

def avgCore(pData,divSet,tcoreList,xcoreList,ycoreList):
    tcore = {}
    xcore = {}
    ycore = {}
    for i in pData.II:
        tcore[i] = np.mean([tcoreList[ll][i] for ll in range(len(tcoreList))])
        for j in pData.Ji[i]:
            xcore[i,j] = np.mean([xcoreList[ll][i,j] for ll in range(len(xcoreList))])
        for par in range(len(divSet[i])):
            ycore[i,par] = np.mean([ycoreList[ll][i,par] for ll in range(len(ycoreList))])
    return tcore,xcore,ycore