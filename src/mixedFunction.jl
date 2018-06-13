# functions of solving mixed first stage problem
function solveMasterMixed(pData,disData,ωInfo,cutSet,tTemp,xTemp,GTemp = Dict())
    M = 400;
    mp = Model(solver = GurobiSolver());
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
      G[i in pData.II,ω in Ω; (i,ω) in ωInfo], Bin
    end);
    @constraint(mp, tFixed[i in pData.II], t[i] == tTemp[i]);
    @constraint(mp, xFixed[i in pData.II, j in pData.Ji[i]], x[i,j] == xTemp[i,j]);
    if GTemp != Dict()
        @constraint(mp, GFixed[(i,ω) in ωInfo], G[i,ω] == GTemp[i,ω]);
    end
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    # logic constraints between G and t
    @constraint(mp, tGcons1[i in pData.II, ω in Ω; (i,ω) in ωInfo], t[i] <= disData[ω].H + G[i,ω]*M);
    @constraint(mp, tGcons2[i in pData.II, ω in Ω; (i,ω) in ωInfo], t[i] >= disData[ω].H - (1 - G[i,ω])*M);

    for ω in Ω
        if cutSet[ω] != []
            for nc in 1:length(cutSet[ω])
                @constraint(mp, θ[ω] >= cutSet[ω][nc][4] + sum(cutSet[ω][nc][1][i]*(mp[:t][i] - cutSet[ω][nc][5][i]) for i in pData.II) +
                    sum(sum(cutSet[ω][nc][2][i,j]*(mp[:x][i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II) +
                    sum(cutSet[ω][nc][3][i,ω1]*(mp[:G][i,ω1] - cutSet[ω][nc][7][i,ω1]) for (i,ω1) in keys(cutSet[ω][nc][7])));
            end
        end
    end
    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));
    return mp;
end
