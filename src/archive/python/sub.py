#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 17:49:36 2019

@author: haoxiangyang
"""

from gurobipy import *

# construct the subproblem for PERT, with y relaxation
def sub_div(pData,dDomega,that,xhat,yhat,divSet,H,M,returnOpt = 0):
    # the subproblem of the divsion algorithm
    # Magnanti-Wong with a small perturbation
    smp = Model('sub_divPrime')
    smp.Params.OutputFlag = 0
    JJ = [(i,j) for i in pData.II for j in pData.Ji[i]]
    x = smp.addVars(JJ, lb = 0, ub = 1, name = "x")
    t = smp.addVars(pData.II, lb = 0, name = "t")
    # relax the logic binary variables
    G = smp.addVars(pData.II, lb = 0, ub = 1, name = "G")
    parI = [(i,par) for i in pData.II for par in range(len(divSet[i]))]
    Gy = smp.addVars(parI, lb = 0, ub = 1, name = "Gy")
    s = smp.addVars(JJ, lb = 0, ub = 1, name = "s")
    
    # add the basic sub problem constraints
    GyRelax1 = smp.addConstrs((Gy[i,par] <= G[i] for i in pData.II for par in range(len(divSet[i]))), name = "GyRelax1")
    GyRelax2 = smp.addConstrs((Gy[i,par] <= yhat[i,par] for i in pData.II for par in range(len(divSet[i]))), name = "GyRelax2")
    GyRelax3 = smp.addConstrs((Gy[i,par] >= G[i] + yhat[i,par] - 1 for i in pData.II for par in range(len((divSet[i])))), name = "GyRelax3")

    GCons1 = smp.addConstrs((dDomega.H*G[i] - sum(H[divSet[i][par].endH]*Gy[i,par] for par in range(len(divSet[i]))) <= dDomega.H - that[i] for i in pData.II), name = "GCons1")
    GCons2 = smp.addConstrs((sum(H[divSet[i][par].startH]*Gy[i,par] for par in range(len(divSet[i]))) - dDomega.H*G[i] >= -that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] \
                             for par in range(len(divSet[i]))) for i in pData.II), name = "GCons2")
    GFixed0 = smp.addConstrs((G[i] >= sum(yhat[i,par] for par in range(len(divSet[i])) if dDomega.H <= H[divSet[i][par].startH]) for i in pData.II), name = "GFixed0")
    GFixed1 = smp.addConstrs((G[i] <= 1 - sum(yhat[i,par] for par in range(len(divSet[i])) if dDomega.H >= H[divSet[i][par].endH]) for i in pData.II), name = "GFixed1")

    # add the predecessors and the successors logic constraints
    GSuccessors = smp.addConstrs((G[i] <= G[k] for i in pData.II for k in pData.Succ[i]), name = "GSuccessors")

    # add the basic sub problem constraints
    tGbound = smp.addConstrs((t[i] >= dDomega.H*G[i] for i in pData.II), name = "tGbound")
    tFnAnt1 = smp.addConstrs((t[i] + G[i]*M[i] >= that[i] for i in pData.II), name = "tFnAnt1")
    tFnAnt2 = smp.addConstrs((t[i] - G[i]*M[i] <= that[i] for i in pData.II), name = "tFnAnt2")
    xFnAnt1 = smp.addConstrs((x[i,j] + G[i] >= xhat[i,j] for i in pData.II for j in pData.Ji[i]), name = "xFnAnt1")
    xFnAnt2 = smp.addConstrs((x[i,j] - G[i] <= xhat[i,j] for i in pData.II for j in pData.Ji[i]), name = "xFnAnt2")

    # linearize the bilinear term of x[i,j]*G[i]
    xGlin1 = smp.addConstrs((s[i,j] <= G[i] for i in pData.II for j in pData.Ji[i]), name = "xGlin1")
    xGlin2 = smp.addConstrs((s[i,j] <= x[i,j] for i in pData.II for j in pData.Ji[i]), name = "xGlin2")
    xGlin3 = smp.addConstrs((s[i,j] >= x[i,j] - 1 + G[i] for i in pData.II for j in pData.Ji[i]), name = "xGlin3")

    budgetConstr = smp.addConstr((sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B), name = "budgetConstr")
    xConstr = smp.addConstrs((sum(x[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II), name = "xConstr")
    durationConstr = smp.addConstrs((t[k[1]] - t[k[0]] >= pData.D[k[0]] + dDomega.d[k[0]]*G[k[0]]
        - sum(pData.D[k[0]]*pData.eff[k[0]][j]*x[k[0],j] + dDomega.d[k[0]]*pData.eff[k[0]][j]*s[k[0],j] for j in pData.Ji[k[0]]) for k in pData.K), name = "durationConstr")

    smp.setObjective(t[0], GRB.MINIMIZE)
    smp.update()
    smp.optimize()
    if smp.Status == 2:
        if returnOpt == 0:
            return smp.objVal
        else:
            return smp.objVal,smp
    else:
        print("Error")

# construct the subproblem for PERT
def sub_divDual(pData,dDomega,that,xhat,yhat,divSet,H,M,tcore,xcore,ycore,returnOpt = 0):
    # the subproblem of the divsion algorithm
    # Magnanti-Wong with a small perturbation
    smp = Model('sub_divPrime')
    smp.Params.OutputFlag = 0
    JJ = [(i,j) for i in pData.II for j in pData.Ji[i]]
    x = smp.addVars(JJ, lb = 0, ub = 1, name = "x")
    t = smp.addVars(pData.II, lb = 0, name = "t")
    # relax the logic binary variables
    G = smp.addVars(pData.II, lb = 0, ub = 1, name = "G")
    parI = [(i,par) for i in pData.II for par in range(len(divSet[i]))]
    Gy = smp.addVars(parI, lb = 0, ub = 1, name = "Gy")
    s = smp.addVars(JJ, lb = 0, ub = 1, name = "s")
    
    # add the basic sub problem constraints
    GyRelax1 = smp.addConstrs((Gy[i,par] <= G[i] for i in pData.II for par in range(len(divSet[i]))), name = "GyRelax1")
    GyRelax2 = smp.addConstrs((Gy[i,par] <= yhat[i,par] for i in pData.II for par in range(len(divSet[i]))), name = "GyRelax2")
    GyRelax3 = smp.addConstrs((Gy[i,par] >= G[i] + yhat[i,par] - 1 for i in pData.II for par in range(len((divSet[i])))), name = "GyRelax3")

    GCons1 = smp.addConstrs((dDomega.H*G[i] - sum(H[divSet[i][par].endH]*Gy[i,par] for par in range(len(divSet[i]))) <= dDomega.H - that[i] for i in pData.II), name = "GCons1")
    GCons2 = smp.addConstrs((sum(H[divSet[i][par].startH]*Gy[i,par] for par in range(len(divSet[i]))) - dDomega.H*G[i] >= -that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] \
                             for par in range(len(divSet[i]))) for i in pData.II), name = "GCons2")
    GFixed0 = smp.addConstrs((G[i] >= sum(yhat[i,par] for par in range(len(divSet[i])) if dDomega.H <= H[divSet[i][par].startH]) for i in pData.II), name = "GFixed0")
    GFixed1 = smp.addConstrs((G[i] <= 1 - sum(yhat[i,par] for par in range(len(divSet[i])) if dDomega.H >= H[divSet[i][par].endH]) for i in pData.II), name = "GFixed1")

    # add the predecessors and the successors logic constraints
    GSuccessors = smp.addConstrs((G[i] <= G[k] for i in pData.II for k in pData.Succ[i]), name = "GSuccessors")

    # add the basic sub problem constraints
    tGbound = smp.addConstrs((t[i] >= dDomega.H*G[i] for i in pData.II), name = "tGbound")
    tFnAnt1 = smp.addConstrs((t[i] + G[i]*M[i] >= that[i] for i in pData.II), name = "tFnAnt1")
    tFnAnt2 = smp.addConstrs((t[i] - G[i]*M[i] <= that[i] for i in pData.II), name = "tFnAnt2")
    xFnAnt1 = smp.addConstrs((x[i,j] + G[i] >= xhat[i,j] for i in pData.II for j in pData.Ji[i]), name = "xFnAnt1")
    xFnAnt2 = smp.addConstrs((x[i,j] - G[i] <= xhat[i,j] for i in pData.II for j in pData.Ji[i]), name = "xFnAnt2")

    # linearize the bilinear term of x[i,j]*G[i]
    xGlin1 = smp.addConstrs((s[i,j] <= G[i] for i in pData.II for j in pData.Ji[i]), name = "xGlin1")
    xGlin2 = smp.addConstrs((s[i,j] <= x[i,j] for i in pData.II for j in pData.Ji[i]), name = "xGlin2")
    xGlin3 = smp.addConstrs((s[i,j] >= x[i,j] - 1 + G[i] for i in pData.II for j in pData.Ji[i]), name = "xGlin3")

    budgetConstr = smp.addConstr((sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B), name = "budgetConstr")
    xConstr = smp.addConstrs((sum(x[i,j] for j in pData.Ji[i]) <= 1 for i in pData.II), name = "xConstr")
    durationConstr = smp.addConstrs((t[k[1]] - t[k[0]] >= pData.D[k[0]] + dDomega.d[k[0]]*G[k[0]]
        - sum(pData.D[k[0]]*pData.eff[k[0]][j]*x[k[0],j] + dDomega.d[k[0]]*pData.eff[k[0]][j]*s[k[0],j] for j in pData.Ji[k[0]]) for k in pData.K), name = "durationConstr")

    smp.setObjective(t[0], GRB.MINIMIZE)
    smp.update()
    smp.optimize()
    vhat = smp.objVal
    Ghat = {}
    for i in pData.II:
        Ghat[i] = G[i].X
        
    # solve the subproblem by dual formulation
    sp = Model('sub_divDual')
    sp.Params.OutputFlag = 0
    λFG1 = sp.addVars(parI, lb = -1e100, ub = 0, name = "λFG1")
    λFG2 = sp.addVars(parI, lb = -1e100, ub = 0, name = "λFG2")
    λFG3 = sp.addVars(parI, lb = 0, name = "λFG3")
    λHG1 = sp.addVars(pData.II, lb = -1e100, ub = 0, name = "λHG1")
    λHG2 = sp.addVars(pData.II, lb = 0, name = "λHG1")
    λGy1 = sp.addVars(pData.II, lb = 0, name = "λGy1")
    λGy2 = sp.addVars(pData.II, lb = -1e100, ub = 0, name = "λGy2")
    λGG = sp.addVars(pData.K, lb = -1e100, ub = 0, name = "λGG")
    λtG1 = sp.addVars(pData.II, lb = 0, name = "λtG1")
    λtG2 = sp.addVars(pData.II, lb = 0, name = "λtG2")
    λtG3 = sp.addVars(pData.II, lb = -1e100, ub = 0, name = "λtG3")
    λxG1 = sp.addVars(JJ, lb = 0, name = "λxG1")
    λxG2 = sp.addVars(JJ, lb = -1e100, ub = 0, name = "λxG2")
    λsG1 = sp.addVars(JJ, lb = -1e100, ub = 0, name = "λsG1")
    λsG2 = sp.addVars(JJ, lb = -1e100, ub = 0, name = "λsG2")
    λsG3 = sp.addVars(JJ, lb = 0, name = "λsG3")
    λbudget = sp.addVar(lb = -1e100, ub = 0, name = "λbudget")
    λxub = sp.addVars(pData.II, lb = -1e100, ub = 0, name = "λxub")
    λdur = sp.addVars(pData.K, lb = 0, name = "λdur")

    tBool = {}
    for i in pData.II:
        if i != 0:
            tBool[i] = 0
        else:
            tBool[i] = 1

    # t constraint
    tConstr = sp.addConstrs((tBool[i] - λtG1[i] - λtG2[i] - λtG3[i] +
        sum(λdur[k] for k in pData.K if k[0] == i) - sum(λdur[k] for k in pData.K if k[1] == i) >= 0 for i in pData.II), name = "tConstr")
    # x constraint
    xConstr = sp.addConstrs((-λxG1[i,j] - λxG2[i,j] + λsG2[i,j] + λsG3[i,j] - pData.b[i][j]*λbudget - \
                             λxub[i] - sum(pData.D[i]*pData.eff[i][j]*λdur[k] for k in pData.K if k[0] == i) >= 0 \
                             for i in pData.II for j in pData.Ji[i]), name = "xConstr")
    # G constraint
    gConstr = sp.addConstrs((sum(λFG1[i,par] + λFG3[i,par] for par in range(len((divSet[i])))) - dDomega.H*(λHG1[i] - λHG2[i]) -
        λGy1[i] - λGy2[i] + sum(λGG[k] for k in pData.K if k[1] == i) - sum(λGG[k] for k in pData.K if k[0] == i) + λtG1[i]*dDomega.H -
        λtG2[i]*M[i] + λtG3[i]*M[i] + sum(-λxG1[i,j] + λxG2[i,j] + λsG1[i,j] + λsG3[i,j] for j in pData.Ji[i]) +
        sum(dDomega.d[i]*λdur[k] for k in pData.K if k[0] == i) >= 0 for i in pData.II), name = "gConstr")
    # Gy constraint
    FConstr = sp.addConstrs((-λFG1[i,par] - λFG2[i,par] - λFG3[i,par] + H[divSet[i][par].endH]*λHG1[i] - \
                            H[divSet[i][par].startH]*λHG2[i] >= 0 for i in pData.II for par in range(len(divSet[i]))), name = "FConstr")
    # s constraint
    sConstr = sp.addConstrs((-λsG1[i,j] - λsG2[i,j] - λsG3[i,j] - sum(dDomega.d[i]*pData.eff[i][j]*λdur[k] for k in pData.K if k[0] == i) >= 0\
                            for i in pData.II for j in pData.Ji[i]),name = "sConstr")

    corePoint = sum(sum(yhat[i,par]*λFG2[i,par] + (yhat[i,par] - 1)*λFG3[i,par] for par in range(len(divSet[i]))) for i in pData.II) + \
        sum(λHG1[i]*(dDomega.H - that[i]) + λHG2[i]*(-that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in range(len(divSet[i])))) for i in pData.II) +\
        sum(λGy1[i]*(sum(yhat[i,par] for par in range(len(divSet[i])) if dDomega.H <= H[divSet[i][par].startH])) +\
            λGy2[i]*(1 - sum(yhat[i,par] for par in range(len(divSet[i])) if dDomega.H >= H[divSet[i][par].endH])) for i in pData.II) +\
        sum(that[i]*(λtG2[i] + λtG3[i]) + sum(xhat[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -\
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +\
        sum(pData.D[k[0]]*λdur[k] for k in pData.K)

    # objective function of the binary feasible solution should be the same
    binaryTight = sp.addConstr(sum(sum(yhat[i,par]*λFG2[i,par] + (yhat[i,par] - 1)*λFG3[i,par] for par in range(len(divSet[i]))) for i in pData.II) +\
        sum(λHG1[i]*(dDomega.H - that[i]) + λHG2[i]*(-that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in range(len(divSet[i])))) for i in pData.II) +\
        sum(λGy1[i]*(sum(yhat[i,par] for par in range(len(divSet[i])) if dDomega.H <= H[divSet[i][par].startH])) +\
            λGy2[i]*(1 - sum(yhat[i,par] for par in range(len(divSet[i])) if dDomega.H >= H[divSet[i][par].endH])) for i in pData.II) +\
        sum(that[i]*(λtG2[i] + λtG3[i]) + sum(xhat[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -\
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +\
        sum(pData.D[k[0]]*λdur[k] for k in pData.K) >= (1 - 1e-5)*vhat)

    # optimize the fractional solution's objective
    sp.setObjective(sum(sum(ycore[i,par]*λFG2[i,par] + (ycore[i,par] - 1)*λFG3[i,par] for par in range(len(divSet[i]))) for i in pData.II) +\
        sum(λHG1[i]*(dDomega.H - tcore[i]) + λHG2[i]*(-tcore[i] + sum(ycore[i,par]*H[divSet[i][par].startH] for par in range(len(divSet[i])))) for i in pData.II) +\
        sum(λGy1[i]*(sum(ycore[i,par] for par in range(len(divSet[i])) if dDomega.H <= H[divSet[i][par].startH])) +\
            λGy2[i]*(1 - sum(ycore[i,par] for par in range(len(divSet[i])) if dDomega.H >= H[divSet[i][par].endH])) for i in pData.II) +\
        sum(tcore[i]*(λtG2[i] + λtG3[i]) + sum(xcore[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -\
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +\
        sum(pData.D[k[0]]*λdur[k] for k in pData.K), GRB.MAXIMIZE)

    sp.update()
    sp.optimize()

    if sp.Status == 2:
        λdict = {}             # dual for x
        πdict = {}             # dual for t
        γdict = {}             # dual for y
        vhat = corePoint.getValue()
        for i in pData.II:
            πdict[i] = -λHG1[i].X - λHG2[i].X + λtG2[i].X + λtG3[i].X
            for j in pData.Ji[i]:
                λdict[i,j] = λxG1[i,j].X + λxG2[i,j].X
            for par in range(len(divSet[i])):
                γdict[i,par] = H[divSet[i][par].startH]*λHG2[i].X + λFG2[i,par].X + λFG3[i,par].X
                if dDomega.H <= H[divSet[i][par].startH]:
                    γdict[i,par] += λGy1[i].X
                elif dDomega.H >= H[divSet[i][par].endH]:
                    γdict[i,par] -= λGy2[i].X
    else:
        vhat = smp.objVal
        λdict = {}             # dual for x
        πdict = {}             # dual for t
        γdict = {}             # dual for y
        for i in pData.II:
            πdict[i] = -GCons1[i].Pi - GCons2[i].Pi + tFnAnt1[i].Pi + tFnAnt2[i].Pi
            for j in pData.Ji[i]:
                λdict[i,j] = xFnAnt1[i,j].Pi + xFnAnt2[i,j].Pi
            for par in range(len(divSet[i])):
                γdict[i,par] = H[divSet[i][par].startH]*GCons2[i].Pi + GyRelax2[i,par].Pi + GyRelax3[i,par].Pi
                if dDomega.H <= H[divSet[i][par].startH]:
                    γdict[i,par] += GFixed0[i].Pi
                elif dDomega.H >= H[divSet[i][par].endH]:
                    γdict[i,par] -= GFixed1[i].Pi  
    
    if returnOpt == 0:
        return πdict,λdict,γdict,vhat,Ghat
    else:
        return πdict,λdict,γdict,vhat,sp