# script for two-stage expected time expected magnitude disruption program
function expModel(pData,eH,ed,M = 9999999)
    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-9));

    @variable(mp,t0[i in pData.II] >= 0);
    @variable(mp,0 <= x0[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(mp,t[i in pData.II] >= 0);
    @variable(mp,0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(mp,G[i in pData.II], Bin);
    @variable(mp,0 <= s[i in pData.II, j in pData.Ji[i]] <= 1);

    @constraint(mp, FConstr[i in pData.II], eH - (1 - G[i])*M <= t0[i]);
    @constraint(mp, GConstr[i in pData.II], eH + G[i]*M >= t0[i]);
    @constraint(mp, tConstr1[i in pData.II], t[i] + G[i]*M >= t0[i]);
    @constraint(mp, tConstr2[i in pData.II], t[i] - G[i]*M <= t0[i]);
    @constraint(mp, tConstr3[i in pData.II], t[i] >= eH * G[i]);
    @constraint(mp, xConstr1[i in pData.II, j in pData.Ji[i]], x[i,j] + G[i] >= x0[i,j]);
    @constraint(mp, xConstr2[i in pData.II, j in pData.Ji[i]], x[i,j] - G[i] <= x0[i,j]);
    @constraint(mp, durationConstr1[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + ed[k[1]]*G[k[1]]
                      - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + ed[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));
    @constraint(mp, durationConstr2[k in pData.K], t0[k[2]] - t0[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x0[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr0[i in pData.II], sum(x0[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(mp, budgetConstr0, sum(sum(pData.b[i][j]*x0[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, Slinear1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(mp, Slinear2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j] + 1 - G[i]);
    @constraint(mp, Slinear3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @objective(mp,Min,pData.p0*t0[0] + (1 - pData.p0)*t[0]);

    solve(mp);
    texp = Dict();
    xexp = Dict();
    gexp = Dict();
    for i in pData.II
        texp[i] = getvalue(mp[:t0][i]);
        for j in pData.Ji[i]
            xexp[i,j] = getvalue(mp[:x0][i,j]);
        end
        gexp[i] = getvalue(mp[:G][i]);
    end
    fexp = getobjectivevalue(mp);
    return texp,xexp,fexp,gexp,mp;
end
