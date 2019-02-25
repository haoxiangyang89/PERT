#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:10:33 2019

@author: haoxiangyang
"""

# auxiliary functions to tighten the branch-and-cut process

import numpy as np

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