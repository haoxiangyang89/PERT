function relaxPart(pData,disData,Ω,cutSet,partCurrent,partDet,M = 9999999)
    # formulate the linear relaxation of the master problem and solve it to obtain the dual solution
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

    lrm = Model(solver = GurobiSolver(OutputFlag = 0));

    @variable(lrm, t[i in pData.II] >= 0);
    @variable(lrm, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(lrm, G[i in pData.II, ω in Ω]);
    @variable(lrm, 0 <= F[i in pData.II, m in 1:partNo[i]] <= 1);
    @variable(lrm, θ[ω in Ω] >= 0);

    @objective(lrm, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    @constraint(lrm, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(lrm, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(lrm, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    # logic constraints between G and t
    @constraint(lrm, tGcons1[i in pData.II, ω in Ω], t[i] <= disData[ω].H + G[i,ω]*M);
    @constraint(lrm, tGcons2[i in pData.II, ω in Ω], t[i] >= disData[ω].H - (1 - G[i,ω])*M);

    # logic constraints of G
    @constraint(lrm, GijCon[i in pData.II, j in pData.Succ[i],ω in Ω], G[i,ω] <= G[j,ω]);
    @constraint(lrm, Gω12Con[i in pData.II, ω in 1:(length(Ω) - 1)], G[i,ω] >= G[i,ω + 1]);
    @constraint(lrm, Gdet1[i in pData.II,ω in Ω; partDet[i][partRev[i][ω]] == 1], G[i,ω] == 1);
    @constraint(lrm, Gdet0[i in pData.II,ω in Ω; partDet[i][partRev[i][ω]] == -1], G[i,ω] == 0);

    # binding constraints of G and F
    @constraint(lrm, GFCon[i in pData.II, ω in Ω], G[i,ω] == F[i,partRev[i][ω]]);

    # add cuts
    @constraint(lrm, cuts[ω in Ω, nc in 1:length(cutSet[ω])], θ[ω] >= cutSet[ω][nc][4] +
        sum(cutSet[ω][nc][1][i]*(lrm[:t][i] - cutSet[ω][nc][5][i]) for i in pData.II) +
        sum(sum(cutSet[ω][nc][3][i,j]*(lrm[:x][i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II) +
        sum(cutSet[ω][nc][2][i]*(lrm[:G][i,ω] - cutSet[ω][nc][7][i]) for i in pData.II));

    # solve the relaxation problem and obtain the dual solution
    solve(lrm);

    μp = Dict();
    for ω in Ω
        for i in pData.II
            μp[i,ω] = getdual(lrm[:GFCon][i,ω]);
        end
    end
    return μp;
end

function createPar(pData,disData,Ω,cutSet,partCurrent)
    # read in the current partition and the cut information to create new partitions
    return partNew;
end
