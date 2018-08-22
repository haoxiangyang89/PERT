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

    tdet = Dict();
    xdet = Dict();
    for i in pData.II
        tdet[i] = getvalue(mp[:t][i]);
        for j in pData.Ji[i]
            xdet[i,j] = getvalue(mp[:x][i,j]);
        end
    end
    fdet = getobjectivevalue(mp);
    return tdet,xdet,fdet;
end
