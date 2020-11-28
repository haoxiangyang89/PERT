function subInt(pData,dDω,xhat,that)
    # solve the MIP recourse problem
    # sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0,Threads = 1));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);

    # add the basic sub problem constraints
    @constraint(sp, fixT[i in pData.II; that[i] < dDω.H], t[i] == that[i]);
    @constraint(sp, boundT[i in pData.II; that[i] >= dDω.H], t[i] >= dDω.H);
    @constraint(sp, fixX[i in pData.II,j in pData.Ji[i]; that[i] < dDω.H], x[i,j] == xhat[i,j]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K; that[k[1]] < dDω.H], t[k[2]] - t[k[1]] >= pData.D[k[1]]*
        (1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(sp, durationConstr2[k in pData.K; that[k[1]] >= dDω.H], t[k[2]] - t[k[1]] >= (pData.D[k[1]]+dDω.d[k[1]])*
        (1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));

    @objective(sp, Min, t[0]);
    solve(sp);
    return getobjectivevalue(sp);
end

function subIntC(pData,dDω,xhat,that,M = 999999,returnOpt = 0)
    # solve the MIP recourse problem
    # sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    global GUROBI_ENV;
    sp = Model(solver = GurobiSolver(GUROBI_ENV,OutputFlag = 0,Threads = 1));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    @variable(sp, G[i in pData.II], Bin);

    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II],dDω.H - (1 - G[i])*M <= that[i]);
    @constraint(sp, GCons[i in pData.II],dDω.H + G[i]*M >= that[i]);

    @constraint(sp, tGcons1[i in pData.II], t[i] <= that[i] + M*G[i]);
    @constraint(sp, tGcons2[i in pData.II], t[i] >= that[i] - M*G[i]);
    @constraint(sp, boundT[i in pData.II], t[i] >= dDω.H*G[i]);
    @constraint(sp, xGcons1[i in pData.II,j in pData.Ji[i]], x[i,j] <= xhat[i,j] + G[i]);
    @constraint(sp, xGcons2[i in pData.II,j in pData.Ji[i]], x[i,j] >= xhat[i,j] - G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(sp, Min, t[0]);
    if returnOpt == 0
        solve(sp);
        return getobjectivevalue(sp);
    else
        return sp;
    end
end

function subIntC_after(pData,dDω,xhat,that,M = 999999,returnOpt = 0)
    # solve the MIP recourse problem
    # sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    global GUROBI_ENV;
    sp = Model(solver = GurobiSolver(GUROBI_ENV,OutputFlag = 0,Threads = 1));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    @variable(sp, G[i in pData.II], Bin);
    @variable(sp, tN >= 0);

    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II],dDω.H - (1 - G[i])*M <= that[i]);
    @constraint(sp, GCons[i in pData.II],dDω.H + G[i]*M >= that[i]);

    @constraint(sp, tGcons1[i in pData.II], t[i] <= that[i] + M*G[i]);
    @constraint(sp, tGcons2[i in pData.II], t[i] >= that[i] - M*G[i]);
    @constraint(sp, boundT[i in pData.II], t[i] >= dDω.H*G[i]);
    @constraint(sp, xGcons1[i in pData.II,j in pData.Ji[i]], x[i,j] <= xhat[i,j] + G[i]);
    @constraint(sp, xGcons2[i in pData.II,j in pData.Ji[i]], x[i,j] >= xhat[i,j] - G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[2]] + dDω.d[k[2]]*G[k[2]]
        - sum(pData.D[k[2]]*pData.eff[k[2]][j]*x[k[2],j] + dDω.d[k[2]]*pData.eff[k[2]][j]*s[k[2],j] for j in pData.Ji[k[2]]));
    @constraint(sp, tNConstr[i in pData.II], tN >= t[i]);

    @objective(sp, Min, tN);
    if returnOpt == 0
        solve(sp);
        return getobjectivevalue(sp);
    else
        return sp;
    end
end

function subIntG(pData,dDω,xhat,that,Ghatω)
    # solve the MIP recourse problem
    M = sum(max(pData.D[i],pData.D[i]+dDω.d[i]) for i in pData.II if i != 0);

    sp = Model(solver = GurobiSolver(OutputFlag = 0,Threads = 1));
    # sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);

    # add the basic sub problem constraints
    @constraint(sp, fixT[i in pData.II; Ghatω[i] == 0], t[i] == that[i]);
    @constraint(sp, boundT[i in pData.II; Ghatω[i] == 1], t[i] >= dDω.H);
    @constraint(sp, fixX[i in pData.II,j in pData.Ji[i]; Ghatω[i] == 0], x[i,j] == xhat[i,j]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K; Ghatω[k[1]] == 0], t[k[2]] - t[k[1]] >= pData.D[k[1]]*
        (1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(sp, durationConstr2[k in pData.K; Ghatω[k[1]] == 1], t[k[2]] - t[k[1]] >= (pData.D[k[1]]+dDω.d[k[1]])*
        (1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));

    @objective(sp, Min, t[0]);
    solve(sp);
    return getobjectivevalue(sp);
end

function subIntGMixed(pData,dDω,xhat,that,ωCurr,Ghatω,M = 400,returnOpt = 0)
    # solve the MIP recourse problem
    #sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0,Threads = 1));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    @variable(sp, s[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, G[i in pData.II], Bin);

    # add the basic sub problem constraints
    @constraint(sp, tGcons1[i in pData.II],dDω.H + G[i]*M >= that[i]);
    @constraint(sp, tGcons2[i in pData.II],dDω.H - (1 - G[i])*M <= that[i]);
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j] + 1 - G[i]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    # with the given Ghat
    @constraint(sp, GConstr[i in pData.II, j in pData.Succ[i]], G[i] <= G[j]);
    @constraint(sp, GGcons1[(i,ω) in keys(Ghatω);ω == ωCurr], G[i] == Ghatω[i,ω]);
    @constraint(sp, GGcons2[(i,ω) in keys(Ghatω);ω < ωCurr], G[i] <= Ghatω[i,ω]);
    @constraint(sp, GGcons3[(i,ω) in keys(Ghatω);ω > ωCurr], G[i] >= Ghatω[i,ω]);

    @objective(sp, Min, t[0]);

    # obtain the dual variables from solving the sub
    if returnOpt == 0
        solve(sp);
        return getobjectivevalue(sp);
    else
        return sp;
    end
end

function ubCal(pData,disData,Ω,xhat,that,bigM,returnOpt = 0)
    # calculate the upper bound of the problem given the master solution
    ubCost = that[0]*pData.p0;
    cω = Dict();
    for ω in Ω
        cω[ω] = subIntC(pData,disData[ω],xhat,that,bigM);
        ubCost += disData[ω].prDis*cω[ω];
    end
    if returnOpt == 0
        return ubCost;
    else
        return ubCost,cω;
    end
end

function ubCalP(pData,disData,Ω,xhat,that,bigM,returnOpt = 0,wp = CachingPool(workers()))
    # parallel version of calculating the upper bound
    ubCost = that[0]*pData.p0;
    cωList = pmap(ω -> subIntC(pData,disData[ω],xhat,that,bigM), wp, Ω);
    ubCost += sum(cωList[i]*disData[Ω[i]].prDis for i in 1:length(Ω));
    if returnOpt == 0
        return ubCost;
    else
        return ubCost,cωList;
    end
end

function ubCalP_after(pData,disData,Ω,xhat,that,bigM,returnOpt = 0,wp = CachingPool(workers()))
    # parallel version of calculating the upper bound
    ubCost = maximum(values(that))*pData.p0;
    cωList = pmap(ω -> subIntC_after(pData,disData[ω],xhat,that,bigM), wp, Ω);
    ubCost += sum(cωList[i]*disData[Ω[i]].prDis for i in 1:length(Ω));
    if returnOpt == 0
        return ubCost;
    else
        return ubCost,cωList;
    end
end

function t_convert(pData,that,xhat,mode)
    # convert the original solution to after solution
    that_new = deepcopy(that);
    xhat_new = deepcopy(xhat);
    if mode == 1
        that_new[0] = 0;
        for i in pData.II
            if i != 0
                that_new[i] = that[i] + pData.D[i]*(1 - sum(pData.eff[i][j]*xhat_new[i,j] for j in pData.Ji[i]));
            end
        end
    elseif mode == 2
        that_new[0] = maximum(values(that_new));
        for i in pData.II
            if i != 0
                that_new[i] = that[i] - pData.D[i]*(1 - sum(pData.eff[i][j]*xhat_new[i,j] for j in pData.Ji[i]));
            end
        end
    end
    return that_new, xhat_new;
end

function ubCalG(pData,disData,Ω,xhat,that,Ghat)
    # calculate the upper bound of the problem given the master solution
    ubCost = that[0]*pData.p0;
    for ω in Ω
        cω = subIntG(pData,disData[ω],xhat,that,Ghat[ω]);
        ubCost += disData[ω].prDis*cω;
    end
    return ubCost;
end
