function createMaster(pData,disData,Ω)
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
    end);
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @objective(mp,Min,pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    return mp;
end

# create a master problem given a partition
function createMaster_Par(pData,disData,Ω,partCurrent,partDet,cutSet = [],M = 9999999)
    # find the matching between scenario and partition
    partRev = Dict();
    partNo = Dict();
    for i in pData.II
        partRev[i] = Dict();
        partNo[i] = length(partCurrent[i]);
        for partIter in 1:partNo[i]
            for item in partCurrent[i][partIter]
                partRev[i][item] = partIter;
            end
        end
    end

    # partCurrent is a list of scenario arrays, each element is a partition of scenarios
    # cutSet includes the cuts generated
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
      G[i in pData.II,m in 1:partNo[i]], Bin
    end);
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    # logic constraints between G and t
    @constraint(mp, tGcons1[i in pData.II, ω in Ω], t[i] <= disData[ω].H + G[i,partRev[i][ω]]*M);
    @constraint(mp, tGcons2[i in pData.II, ω in Ω], t[i] >= disData[ω].H - (1 - G[i,partRev[i][ω]])*M);

    # logic constraints of G
    @constraint(mp, GijCon[i in pData.II, j in pData.Succ[i],ω in Ω], G[i,partRev[i][ω]] <= G[j,partRev[i][ω]]);
    @constraint(mp, Gω12Con[i in pData.II, m in 1:(partNo[i] - 1)], G[i,m] <= G[i,m + 1]);
    @constraint(mp, Gdet1[i in pData.II,m in 1:partNo[i]; partDet[i][m] == 1], G[i,m] == 1);
    @constraint(mp, Gdet0[i in pData.II,m in 1:partNo[i]; partDet[i][m] == -1], G[i,m] == 0);
    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    return mp;
end
