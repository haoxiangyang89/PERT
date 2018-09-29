function createMaster(pData,disData,Ω,Tmax = 999999)
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      0 <= t[i in pData.II] <= Tmax
    end);
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @objective(mp,Min,pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    return mp;
end

# create a master problem given a partition
function createMaster_Par(pData,disData,Ω,partCurrent,partDet,cutSet,M = 9999999)
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
    mp = Model(solver = GurobiSolver());
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
    @constraint(mp, GijCon[i in pData.II, j in pData.Succ[i],ω in Ω], G[i,partRev[i][ω]] <= G[j,partRev[j][ω]]);
    @constraint(mp, Gω12Con[i in pData.II, m in 1:(partNo[i] - 1)], G[i,m] >= G[i,m + 1]);
    @constraint(mp, Gdet1[i in pData.II,m in 1:partNo[i]; partDet[i][m] == 1], G[i,m] == 1);
    #@constraint(mp, Gdet0[i in pData.II,m in 1:partNo[i]; partDet[i][m] == -1], G[i,m] == 0);

    for ω in Ω
        if cutSet[ω] != []
            for nc in 1:length(cutSet[ω])
                @constraint(mp, θ[ω] >= cutSet[ω][nc][4] + sum(cutSet[ω][nc][1][i]*(mp[:t][i] - cutSet[ω][nc][5][i]) for i in pData.II) +
                    sum(sum(cutSet[ω][nc][3][i,j]*(mp[:x][i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II) +
                    sum(cutSet[ω][nc][2][i]*(mp[:G][i,partRev[i][ω]] - cutSet[ω][nc][7][i]) for i in pData.II));
            end
        end
    end
    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    return mp;
end

# create a master problem with a selection of binary variables
function createMaster_Mixed(pData,disData,Ω,ωInfo,cutSet,M = 9999999)
    H = Dict();
    H[0] = 0;
    H[length(Ω)+1] = M;
    for ω in Ω
        H[ω] = disData[ω].H;
    end
    # cutSet includes the cuts generated
    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-9,FeasibilityTol = 1e-9));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
      G[i in pData.II,ω in Ω; (i,ω) in ωInfo], Bin
    end);
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

# create a master problem with a selection of binary variables
function createMaster_MixedL(pData,disData,Ω,ωInfo,cutSet,partCurrentTemp,partDetTemp,M = 999999)
    H = Dict();
    H[0] = 0;
    H[length(Ω)+1] = M;
    for ω in Ω
        H[ω] = disData[ω].H;
    end
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
    # cutSet includes the cuts generated
    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-9,FeasibilityTol = 1e-9));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
      0 <= G[i in pData.II,ω in Ω] <= 1
    end);
    for (i,ω) in ωInfo
        setcategory(G[i,ω], :Bin);
    end
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    # logic constraints between G and t
    @constraint(mp, tGcons1I[i in pData.II, ω in Ω], t[i] <= H[ω] + G[i,ω]*M);
    @constraint(mp, tGcons2I[i in pData.II, ω in Ω], t[i] >= H[ω] - (1 - G[i,ω])*M);
    @constraint(mp, tGcons3I[i in pData.II], t[i] >= sum(G[i,ω]*(H[ω] - H[ω - 1]) for ω in Ω));
    # logic constraints of G
    @constraint(mp, GijCon[i in pData.II, j in pData.Succ[i],ω in Ω], G[i,ω] <= G[j,ω]);
    @constraint(mp, Gω12Con[i in pData.II, ω in 1:(length(Ω)-1)], G[i,ω] >= G[i,ω + 1]);
    @constraint(mp, Gdet1[i in pData.II, ω in Ω; partDet[i][partRev[i][ω]] == 1], G[i,ω] == 1);
    @constraint(mp, Gdet0[i in pData.II, ω in Ω; partDet[i][partRev[i][ω]] == -1], G[i,ω] == 0);

    for ω in Ω
        if cutSet[ω] != []
            for nc in 1:length(cutSet[ω])
                @constraint(mp, θ[ω] >= cutSet[ω][nc][4] + sum(cutSet[ω][nc][1][i]*(mp[:t][i] - cutSet[ω][nc][5][i]) for i in pData.II) +
                    sum(sum(cutSet[ω][nc][2][i,j]*(mp[:x][i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II) +
                    sum(cutSet[ω][nc][3][i]*(mp[:G][i,ω] - cutSet[ω][nc][7][i]) for i in pData.II));
            end
        end
    end

    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    return mp;
end

function masterF(pData,disData,Ω,ωInfo,cutSet,Tmax = 999999)
    H = Dict();
    H[0] = 0;
    H[length(Ω)+1] = Tmax;
    for ω in Ω
        H[ω] = disData[ω].H;
    end

    # the master with generated columns
    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-9,FeasibilityTol = 1e-9));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
      0 <= F[i in pData.II,ω in 0:length(Ω)] <= 1
    end);
    for (i,ω) in ωInfo
        setcategory(F[i,ω], :Bin);
    end

    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    # logic constraints between G and t
    @constraint(mp, tFcons1I[i in pData.II], t[i] <= sum(H[ω+1]*F[i,ω] for ω in Ω));
    @constraint(mp, tFcons2I[i in pData.II], t[i] >= sum(H[ω]*F[i,ω] for ω in Ω));
    # logic constraints of G
    @constraint(mp, Fsum[i in pData.II], sum(F[i,ω] for ω in 0:length(Ω)) == 1);

    for ω in Ω
        if cutSet[ω] != []
            for nc in 1:length(cutSet[ω])
                @constraint(mp, θ[ω] >= cutSet[ω][nc][4] + sum(cutSet[ω][nc][1][i]*(mp[:t][i] - cutSet[ω][nc][5][i]) for i in pData.II) +
                    sum(sum(cutSet[ω][nc][2][i,j]*(mp[:x][i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II) +
                    sum(sum(cutSet[ω][nc][3][i]*(mp[:F][i,ω1] - cutSet[ω][nc][7][i,ω1]) for ω1 in Ω if ω1 >= ω) for i in pData.II));
            end
        end
    end

    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));
    return mp;
end

function updateMaster(mp,ubInfo,lbInfo)
    for i in pData.II
        @constraint(mp,mp[:t][i] <= ubInfo[i]);
        @constraint(mp,mp[:t][i] >= lbInfo[i]);
    end
    return mp;
end

function createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax,yCuts = 0)
    H = Dict();
    H[0] = 0;
    H[length(Ω)+1] = Tmax;
    for ω in Ω
        H[ω] = disData[ω].H;
    end

    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-9,FeasibilityTol = 1e-9,OutputFlag = 0));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
      y[i in pData.II, par in 1:length(divSet[i])], Bin
    end);
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(mp, tub[i in pData.II], t[i] <= sum(H[divSet[i][par].endH]*y[i,par] for par in 1:length(divSet[i])));
    @constraint(mp, tlb[i in pData.II], t[i] >= sum(H[divSet[i][par].startH]*y[i,par] for par in 1:length(divSet[i])));
    @constraint(mp, yConstr[i in pData.II], sum(y[i,par] for par in 1:length(divSet[i])) == 1);
    @constraint(mp, yLimit[i in pData.II, par in 1:length(divSet[i]); divDet[i][par] != 0], y[i,par] == 0);
    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    # add the cut
    # cutInfo = 2 dimensional vector, first dimention record the primal solution,
    # second dimension record the dual solution for every scenario
    for nc in 1:length(cutSet)
        for ω in Ω
            if cutSet[nc][2][ω] != []
                # if there is a cut for scenario ω in iteration nc
                # @constraint(mp, θ[ω] >= cutSet[nc][2][ω][1] + sum(cutSet[nc][2][ω][2][i]*(mp[:t][i] - cutSet[nc][1][1][i]) for i in pData.II) +
                #     sum(sum(cutSet[nc][2][ω][3][i,j]*(mp[:x][i,j] - cutSet[nc][1][2][i,j]) for j in pData.Ji[i]) for i in pData.II) +
                #     sum(sum(cutSet[nc][2][ω][4][i,revPar(cutSet[nc][1][4][i],divSet[i][par])]*(mp[:y][i,par] - cutSet[nc][1][3][i,revPar(cutSet[nc][1][4][i],divSet[i][par])]) for par in 1:length(divSet[i])) for i in pData.II));
                @constraint(mp, θ[ω] >= cutSet[nc][2][ω][1] + sum(cutSet[nc][2][ω][2][i]*(mp[:t][i] - cutSet[nc][1][1][i]) for i in pData.II) +
                    sum(sum(cutSet[nc][2][ω][3][i,j]*(mp[:x][i,j] - cutSet[nc][1][2][i,j]) for j in pData.Ji[i]) for i in pData.II) +
                    sum(sum(cutSet[nc][2][ω][4][i,par]*(sum(mp[:y][i,parNew] for parNew in 1:length(divSet[i]) if revPar(cutSet[nc][1][4][i],divSet[i][parNew]) == par) - cutSet[nc][1][3][i,par]) for par in 1:length(cutSet[nc][1][4][i])) for i in pData.II));
            end
        end
    end

    # add the constraints between y
    if yCuts != 0
        for k in pData.K
            # for each precedence relationship
            for par1 in 1:length(divSet[k[1]])
                for par2 in 1:length(divSet[k[2]])
                    if H[divSet[k[2]][par2].endH] < H[divSet[k[1]][par1].startH] + pData.D[k[1]]*(1 - maximum(values(pData.eff[k[1]])))
                        @constraint(mp, y[k[1],par1] + y[k[2],par2] <= 1);
                    end
                end
            end
        end
    end

    return mp;
end
