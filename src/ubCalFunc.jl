function subInt(pData,dDω,xhat,that)
    # solve the MIP recourse problem
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
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

function subIntG(pData,dDω,xhat,that,Ghatω)
    # solve the MIP recourse problem
    M = sum(max(pData.D[i],pData.D[i]+dDω.d[i]) for i in pData.II if i != 0);

    # sp = Model(solver = GurobiSolver(OutputFlag = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
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

function ubCal(pData,disData,Ω,xhat,that)
    # calculate the upper bound of the problem given the master solution
    ubCost = that[0]*pData.p0;
    for ω in Ω
        cω = subInt(pData,disData[ω],xhat,that);
        ubCost += disData[ω].prDis*cω;
    end
    return ubCost;
end

function ubCalP(pData,disData,Ω,xhat,that)
    # parallel version of calculating the upper bound
    ubCost = that[0]*pData.p0;
    cωList = pmap(ω -> subInt(pData,disData[ω],xhat,that), Ω);
    ubCost += sum(cωList[i]*disData[Ω[i]].prDis for i in 1:length(ω));
    return ubCost;
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
