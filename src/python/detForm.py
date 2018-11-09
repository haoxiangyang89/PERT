#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 22:55:59 2018

@author: haoxiangyang
"""

# deterministic formulation
from gurobipy import *

def detBuild(pData):
    # build the deterministic crashing optimization problem
    mp = Model('deterministic')
    mp.Params.OutputFlag = 0

    t = mp.addVars(pData.II)
    xInd = []
    for i in pData.II:
        for j in pData.Ji[i]:
            xInd.append((i,j))
    x = mp.addVars(xInd, ub = 1)

    durationConstr = mp.addConstrs(t[k[1]] - t[k[0]] >= pData.D[k[0]]*(1 - sum(pData.eff[k[0]][j]*x[k[0],j] 
        for j in pData.Ji[k[0]])) for k in pData.K)
    budgetConstr = mp.addConstr(sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B)
    xConstr = mp.addConstrs(sum(x[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II)

    mp.setObjective(t[0], GRB.MINIMIZE)
    mp.optimize()

    tdet = {}
    xdet = {}
    for i in pData.II:
        tdet[i] = t[i].X
        for j in pData.Ji[i]:
            xdet[i,j] = x[i,j].X
    fdet = mp.objVal
    return tdet,xdet,fdet