# build the LP to obtain the earliest possible starting time of an activity
function iSolve(pData,iTarget,ubTemp = 9999999)
    # mp = Model(solver = GurobiSolver(OutputFlag = 0));
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(mp, t[i in pData.II] >= 0);
    @variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);

    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 -
        sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    @constraint(mp, ubCon, t[0] <= ubTemp);

    @objective(mp, Min, t[iTarget]);
    solve(mp);

    return getobjectivevalue(mp);
end

# build the LP to obtain the latest possible starting time of an activity
function lSolve(pData,iTarget,ubTemp = 999999)
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(mp, t[i in pData.II] >= 0);
    @variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);

    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 -
        sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    @constraint(mp, ubCon, t[0] <= ubTemp);

    @objective(mp, Max, t[iTarget]);
    solve(mp);

    return getobjectivevalue(mp);
end

# build the LP to obtain the big M
function iSolve_NC(pData,dDω,iTarget,brInfoω)
    # mp = Model(solver = GurobiSolver(OutputFlag = 0));
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(mp, t[i in pData.II] >= 0);
    Dω = Dict();

    for i in keys(dDω.d)
        if brInfoω[findin(pData.II,i)[1]] == 0
            Dω[i] = max(0,dDω.d[i])+pData.D[i];
        elseif brInfoω[findin(pData.II,i)[1]] == -1
            Dω[i] = pData.D[i];
        else
            Dω[i] = pData.D[i] + dDω.d[i];
        end
    end
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= Dω[k[1]]);

    @objective(mp, Min, t[iTarget]);
    solve(mp);

    return getobjectivevalue(mp);
end

# build the LP to obtain the earliest possible starting time of an activity, after-start disruption
function iSolve_after(pData,iTarget,ubTemp = 9999999)
    # mp = Model(solver = GurobiSolver(OutputFlag = 0));
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(mp, t[i in pData.II] >= 0);
    @variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);

    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[2]]*(1 -
        sum(pData.eff[k[2]][j]*x[k[2],j] for j in pData.Ji[k[2]])));
    @constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    @constraint(mp, ubCon[i in pData.II], t[i] <= ubTemp);

    @objective(mp, Min, t[iTarget]);
    solve(mp);

    return getobjectivevalue(mp);
end

# build the LP to obtain the latest possible starting time of an activity, after-start disruption
function lSolve_after(pData,iTarget,ubTemp = 999999)
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(mp, t[i in pData.II] >= 0);
    @variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);

    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[2]]*(1 -
        sum(pData.eff[k[2]][j]*x[k[2],j] for j in pData.Ji[k[2]])));
    @constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    @constraint(mp, ubCon[i in pData.II], t[i] <= ubTemp);

    @objective(mp, Max, t[iTarget]);
    solve(mp);

    return getobjectivevalue(mp);
end
