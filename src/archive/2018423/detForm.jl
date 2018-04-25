# deterministic/expected counterpart for testing, with crashing options continuous

function detBuild(pData)
    # build the deterministic crashing optimization problem
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    # mp = Model(solver = ClpSolver());

    @variable(mp, t[i in pData.II] >= 0);
    @variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);

    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 -
        sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    @objective(mp, Min, t[0]);
    solve(mp);

    tdet = getvalue(mp[:t]);
    xdet = getvalue(mp[:x]);
    fdet = getobjectivevalue(mp);
    return tdet,xdet,fdet;
end

function expBuild(pData,disData,Ω)
    # build the expected crashing optimization problem
    # the disruption time and the disruption magnitude are both at their expected values

    # obtain the expected value of the disruption time
    expH = sum(disData[ω].H for ω in Ω)/length(Ω);

    # obtain the expected disruption magnitude
    expd = Dict();
    for i in pData.II
        if i != 0
            expd[i] = sum(disData[ω].d[i] for ω in Ω)/length(Ω);
        end
    end
    M = sum(max(pData.D[i],pData.D[i]+expd[i]) for i in pData.II if i != 0);

    # set up the expected crashing optimization problem
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    # mp = Model(solver = ClpSolver());
    @variable(mp, t[i in pData.II] >= 0);
    @variable(mp, t1[i in pData.II] >= 0);
    @variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(mp, 0 <= x1[i in pData.II, j in pData.Ji[i]] <= 1);

    @variable(mp, F[i in pData.II], Bin);
    @variable(mp, G[i in pData.II], Bin);
    @variable(mp, 0 <= S[i in pData.II, j in pData.Ji[i]] <= 1);

    @constraint(mp, FCons[i in pData.II], t[i] >= expH - F[i]*M);
    @constraint(mp, GCons[i in pData.II], t[i] <= expH - 1e-6 + G[i]*M);
    @constraint(mp, FGCons[i in pData.II], F[i] + G[i] == 1);

    @constraint(mp, FtCons1[i in pData.II], t1[i] <= t[i] + M*(1 - F[i]));
    @constraint(mp, FtCons2[i in pData.II], t1[i] >= t[i] - M*(1 - F[i]));
    @constraint(mp, GtCons[i in pData.II], t1[i] >= expH*G[i]);
    @constraint(mp, FxCons1[i in pData.II, j in pData.Ji[i]], x1[i,j] <= x[i,j] + (1 - F[i]));
    @constraint(mp, FxCons2[i in pData.II, j in pData.Ji[i]], x1[i,j] >= x[i,j] - (1 - F[i]));

    @constraint(mp, SCons1[i in pData.II, j in pData.Ji[i]], S[i,j] <= G[i]);
    @constraint(mp, SCons2[i in pData.II, j in pData.Ji[i]], S[i,j] <= x1[i,j] + 1 - G[i]);
    @constraint(mp, SCons3[i in pData.II, j in pData.Ji[i]], S[i,j] <= x1[i,j] - 1 + G[i]);

    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 -
        sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, durationConstr1[k in pData.K], t1[k[2]] - t1[k[1]] >= pData.D[k[1]] + expd[k[1]]*G[k[1]] -
        pData.D[k[1]]*sum(pData.eff[k[1]][j]*x1[k[1],j] for j in pData.Ji[k[1]]) -
        expd[k[1]]*sum(pData.eff[k[1]][j]*S[k[1],j] for j in pData.Ji[k[1]]));
    @constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(mp, budgetConstr1, sum(sum(x1[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr1[i in pData.II], sum(x1[i,j] for j in pData.Ji[i]) <= 1);

    @objective(mp, Min, pData.p0*t[0] + (1 - pData.p0)*t1[0]);
    solve(mp);

    texp = getvalue(mp[:t]);
    xexp = getvalue(mp[:x]);
    fexp = getobjectivevalue(mp);
    return texp,xexp,fexp;
end
