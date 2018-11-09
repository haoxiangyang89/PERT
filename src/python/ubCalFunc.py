#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 01:32:05 2018

@author: haoxiangyang
"""

from gurobipy import *

# functions to calculate the upper bound

def subIntC(pData,dDomega,xhat,that,M,returnOpt = 0):
    # solve the MIP recourse problem
    sp = Model()
    sp.Params.OutputFlag = 0
    t = sp.addVars(pData.II)
    xInd = []
    for i in pData.II:
        for j in pData.Ji[i]:
            xInd.append((i,j))
    x = sp.addVars(xInd, ub = 1)
    s = sp.addVars(xInd, ub = 1)
    G = sp.addVars(pData.II, vtype = GRB.BINARY)

    # add the basic sub problem constraints
    FCons = sp.addConstrs(dDomega.H - (1 - G[i])*M <= that[i] for i in pData.II)
    GCons = sp.addConstrs(dDomega.H + G[i]*M >= that[i] for i in pData.II)

    tGcons1 = sp.addConstrs(t[i] <= that[i] + M*G[i] for i in pData.II)
    tGcons2 = sp.addConstrs(t[i] >= that[i] - M*G[i] for i in pData.II)
    boundT = sp.addConstrs(t[i] >= dDomega.H*G[i] for i in pData.II)
    xGcons1 = sp.addConstrs(x[i,j] <= xhat[i,j] + G[i] for i in pData.II for j in pData.Ji[i])
    xGcons2 = sp.addConstrs(x[i,j] >= xhat[i,j] - G[i] for i in pData.II for j in pData.Ji[i])

    budgetConstr = sp.addConstr(sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B)
    xConstr = sp.addConstrs(sum(x[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II)
    # linearize the bilinear term of x[i,j]*G[i]
    xGlin1 = sp.addConstrs(s[i,j] <= G[i] for i in pData.II for j in pData.Ji[i])
    xGlin2 = sp.addConstrs(s[i,j] <= x[i,j] for i in pData.II for j in pData.Ji[i])
    xGlin3 = sp.addConstrs(s[i,j] >= x[i,j] - 1 + G[i] for i in pData.II for j in pData.Ji[i])

    durationConstr = sp.addConstrs(t[k[1]] - t[k[0]] >= pData.D[k[0]] + dDomega.d[k[0]]*G[k[0]]
        - sum(pData.D[k[0]]*pData.eff[k[0]][j]*x[k[0],j] + dDomega.d[k[0]]*pData.eff[k[0]][j]*s[k[0],j] for j in pData.Ji[k[0]]) for k in pData.K)

    sp.setObjective(t[0], GRB.MINIMIZE)
    if returnOpt == 0:
        sp.optimize()
        return sp.objVal
    else:
        return sp

def ubCal(pData,disData,Omega,xhat,that,bigM,returnOpt = 0):
    ubCost = that[0]*pData.p0
    comega = {}
    for omega in Omega:
        comega[omega] = subIntC(pData,disData[omega],xhat,that,bigM)
        ubCost += disData[omega].prDis*comega[omega]
    if returnOpt == 0:
        return ubCost
    else:
        return ubCost,comega