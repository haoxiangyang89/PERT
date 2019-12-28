#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 11:52:38 2019

@author: haoxiangyang
"""

from multiprocessing import Pool
from gurobipy import *

def subInt(pData,dDomega,xh,th,M = 999999,returnOpt = 0):
    # solve the mixed integer subproblem to compute the scenario ub
    sp = Model('subInt')
    sp.Params.OutputFlag = 0
    t = sp.addVars(pData.II, lb = 0)
    JJ = [(i,j) for i in pData.II for j in pData.Ji[i]]
    x = sp.addVars(JJ, lb = 0, ub = 1)
    s = sp.addVars(JJ, lb = 0, ub = 1)
    G = sp.addVars(pData.II, vtype = GRB.BINARY)

    # add the basic sub problem constraints
    FCons = sp.addConstrs(dDomega.H - (1 - G[i])*M <= th[i] for i in pData.II)
    GCons = sp.addConstrs(dDomega.H + G[i]*M >= th[i] for i in pData.II)

    tGcons1 = sp.addConstrs(t[i] <= th[i] + M*G[i] for i in pData.II)
    tGcons2 = sp.addConstrs(t[i] >= th[i] - M*G[i] for i in pData.II)
    boundT = sp.addConstrs(t[i] >= dDomega.H*G[i] for i in pData.II)
    xGcons1 = sp.addConstrs(x[i,j] <= xh[i,j] + G[i] for i in pData.II for j in pData.Ji[i])
    xGcons2 = sp.addConstrs(x[i,j] >= xh[i,j] - G[i] for i in pData.II for j in pData.Ji[i])

    budgetConstr = sp.addConstr(sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B)
    xConstr = sp.addConstrs(sum(x[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II for j in pData.Ji[i])
    # linearize the bilinear term of x[i,j]*G[i]
    xGlin1 = sp.addConstrs(s[i,j] <= G[i] for i in pData.II for j in pData.Ji[i])
    xGlin2 = sp.addConstrs(s[i,j] <= x[i,j] for i in pData.II for j in pData.Ji[i])
    xGlin3 = sp.addConstrs(s[i,j] >= x[i,j] - 1 + G[i] for i in pData.II for j in pData.Ji[i])

    durationConstr = sp.addConstrs(t[k[1]] - t[k[0]] >= pData.D[k[0]] + dDomega.d[k[0]]*G[k[0]]\
        - sum(pData.D[k[0]]*pData.eff[k[0]][j]*x[k[0],j] + dDomega.d[k[0]]*pData.eff[k[0]][j]*s[k[0],j] for j in pData.Ji[k[0]]) \
        for k in pData.K)

    sp.setObjective(t[0],GRB.MINIMIZE)
    sp.update()
    if returnOpt == 0:
        sp.optimize()
        return sp.objVal
    else:
        return sp
    
def ubCal(pData,disData,Omega,xhat,that,bigM,returnOpt = 0):
    # calculate the upper bound of the problem given the master solution
    ubCost = that[0]*pData.p0
    comega = {}
    for omega in Omega:
        comega[omega] = subInt(pData,disData[omega],xhat,that,bigM)
        ubCost += disData[omega].prDis*comega[omega]
    if returnOpt == 0:
        return ubCost
    else:
        return ubCost,comega

def ubCalP(pData,disData,Omega,xhat,that,bigM,threads = 20,returnOpt = 0):
    # parallel version of calculating the upper bound
    ubCost = that[0]*pData.p0
    p = Pool(threads)
    ans = p.map(functools.partial(subInt,xhat,that), range(50))