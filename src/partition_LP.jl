function relaxPart(pData,disData,Ω,cutSet,partCurrentTemp,partDet,M = 9999999,FPre = [],returnOpt = 0)
    # formulate the linear relaxation of the master problem and solve it to obtain the dual solution
    partRev = Dict();
    partNo = Dict();
    for i in pData.II
        partRev[i] = Dict();
        partNo[i] = length(partCurrentTemp[i]);
        for partIter in 1:partNo[i]
            for item in partCurrentTemp[i][partIter]
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
#    @constraint(lrm, Gdet0[i in pData.II,ω in Ω; partDet[i][partRev[i][ω]] == -1], G[i,ω] == 0);

    # binding constraints of G and F
    if FPre != []
        @constraint(lrm, GFCon[i in pData.II, ω in Ω], G[i,ω] == FPre[i,partRev[i][ω]]);
    else
        @constraint(lrm, GFCon[i in pData.II, ω in Ω], G[i,ω] == F[i,partRev[i][ω]]);
    end

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
            μp[i,ω] = -getdual(lrm[:GFCon][i,ω]);
        end
    end
    if returnOpt == 0
        return μp;
    else
        return μp,lrm;
    end
end

function createPar(pData,disData,Ω,partCurrentTemp,partDetTemp,μp)
    # read in the current partition and the cut information to create new partitions
    partNew = Dict();
    for i in pData.II
        partNew[i] = [];
        for pci in 1:length(partCurrentTemp[i])
            push!(partNew[i],partCurrentTemp[i][pci]);
        end
    end
    for i in pData.II
        pciBest = 1;
        μFPBest = 0;
        μFPBestInd = 0;
        splitBool = false;
        for pci in 1:length(partNew[i])
            μFCurrent = sum(μp[i,ω] for ω in partNew[i][pci]);
            μFBest = μFCurrent;
            μFBestInd = length(partNew[i][pci]);
            # find the best spliting point for the current partition
            if partDet[i][pci] == 0
                for ind in 1:(length(partNew[i][pci]) - 1)
                    firstP = max(0,sum(μp[i,ωi] for ωi in partNew[i][pci][1:ind]));
                    secondP = max(0,sum(μp[i,ωi] for ωi in partNew[i][pci][(ind+1):length(partNew[i][pci])]));
                    if firstP + secondP > μFBest
                        μFBest = firstP + secondP;
                        μFBestInd = ind;
                    end
                end
                if μFBest - μFCurrent > μFPBest
                    μFPBest = μFBest - μFCurrent;
                    μFPBestInd = μFBestInd;
                    pciBest = pci;
                    splitBool = true;
                end
            end
        end
        # split according to the obtained μFPBestInd
        if splitBool
            splice!(partNew[i],pciBest,[partNew[i][pciBest][1:μFPBestInd],partNew[i][pciBest][μFPBestInd + 1:length(partNew[i][pciBest])]]);
            splice!(partDet[i],pciBest,[0,0]);
        end
    end

    return partNew;
end

function createParCut(pData,disData,Ω,partCurrentTemp,partDetTemp,μp)
    # read in the current partition and the cut information to create new partitions
    partNew = Dict();
    partDetNew = Dict();
    for i in pData.II
        partNew[i] = [];
        partDetNew[i] = [];
    end
    for i in pData.II
        for pci in 1:length(partCurrentTemp[i])
            splitBool = false;
            # find the best spliting point for the current partition
            if partDetTemp[i][pci] == 0
                μFPBestInd = 1;
                for ω in 2:length(partCurrentTemp[i][pci])
                    if (μp[i,partCurrentTemp[i][pci][ω-1]] > 0)&&(μp[i,partCurrentTemp[i][pci][ω]] <= 0)
                        splitBool = true;
                        μFPBestInd = ω - 1;
                    end
                end
                # split according to the obtained μFPBestInd
                if splitBool
                    push!(partNew[i],partCurrentTemp[i][pci][1:μFPBestInd])
                    push!(partNew[i],partCurrentTemp[i][pci][μFPBestInd + 1:length(partCurrentTemp[i][pci])]);
                    push!(partDetNew[i],0);
                    push!(partDetNew[i],0);
                else
                    push!(partNew[i],partCurrentTemp[i][pci]);
                    push!(partDetNew[i],partDetTemp[i][pci]);
                end
            else
                push!(partNew[i],partCurrentTemp[i][pci]);
                push!(partDetNew[i],partDetTemp[i][pci]);
            end
        end
    end

    return partNew,partDetNew;
end

function createParOpt(pData,disData,Ω,partCurrentTemp,partDetTemp,cutSet,Tmax = 200)
    # formulate the linear relaxation of the master problem and solve it to obtain the dual solution
    partRev = Dict();
    partNo = Dict();
    for i in pData.II
        partRev[i] = Dict();
        partNo[i] = length(partCurrentTemp[i]);
        for partIter in 1:partNo[i]
            for item in partCurrentTemp[i][partIter]
                partRev[i][item] = partIter;
            end
        end
    end
    # obtain the range where the F is undecided
    partZero = Dict();
    for i in pData.II
        zeroMin = minimum([ω for ω in Ω if partDetTemp[i][partRev[i][ω]] == 0]) - 1;
        zeroMax = maximum([ω for ω in Ω if partDetTemp[i][partRev[i][ω]] == 0]);
        partZero[i] = zeroMin:zeroMax;
    end
    H = Dict();
    H[0] = 0;
    H[length(Ω)+1] = Tmax;
    for ω in Ω
        H[ω] = disData[ω].H;
    end

    mpar = Model(solver = GurobiSolver(IntFeasTol = 1e-9,FeasibilityTol = 1e-9));
    @variable(mpar,t[i in pData.II] >= 0);
    @variable(mpar,0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(mpar,θ[ω in Ω] >= 0);
    @variable(mpar,0 <= G[i in pData.II, ω in Ω] <= 1);
    #@variable(mpar,F[i in pData.II, ω in partZero[i]], Bin);
    @variable(mpar,0 <= F[i in pData.II, ω in partZero[i]] <= 1);

    @constraint(mpar, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mpar, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mpar, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    # logic constraints between G and t
    @constraint(mpar, tGcons1I[i in pData.II], t[i] <= sum(H[ω+1]*F[i,ω] for ω in partZero[i]));
    @constraint(mpar, tGcons2I[i in pData.II], t[i] >= sum(H[ω]*F[i,ω] for ω in partZero[i]));
    @constraint(mpar, Fcons[i in pData.II], sum(F[i,ω] for ω in partZero[i]) == 1);
    @constraint(mpar, GFcons[i in pData.II, ω in Ω], G[i,ω] == sum(F[i,ω1] for ω1 in partZero[i] if ω1 >= ω));
    @constraint(mpar, Gdet1[i in pData.II,ω in Ω; partDetTemp[i][partRev[i][ω]] == 1], G[i,ω] == 1);

    @constraint(mpar, cuts[ω in Ω, nc in 1:length(cutSet[ω])], θ[ω] >= cutSet[ω][nc][4] +
        sum(cutSet[ω][nc][1][i]*(mpar[:t][i] - cutSet[ω][nc][5][i]) for i in pData.II) +
        sum(sum(cutSet[ω][nc][3][i,j]*(mpar[:x][i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II) +
        sum(cutSet[ω][nc][2][i]*(mpar[:G][i,ω] - cutSet[ω][nc][7][i]) for i in pData.II));

    @objective(mpar, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    solve(mpar);
    Fsol = getvalue(mpar[:F]);
    splitLoc = Dict();
    for i in pData.II
        for ω in Ω
            if Fsol[i,ω] == 1
                splitLoc[i] = ω;
            end
        end
    end
    partNew = Dict();
    partDetNew = Dict();
    for i in pData.II
        partNew[i] = [];
        partDetNew[i] = [];
        if i in keys(splitLoc)
            for pci in 1:length(partCurrentTemp[i])
                if !(splitLoc[i] in partCurrentTemp[i][pci])
                    push!(partNew[i],partCurrentTemp[i][pci]);
                    push!(partDetNew[i],partDetTemp[i][pci]);
                else
                    splitPos = findfirst(partCurrentTemp[i][pci],splitLoc[i]);
                    if (splitPos != length(partDetTemp[i][pci]))
                        partNew1 = [partCurrentTemp[i][pci][ps] for ps in 1:splitPos];
                        partNew2 = [partCurrentTemp[i][pci][ps] for ps in (splitPos + 1):length(partCurrentTemp[i][pci])];
                    end
                    push!(partNew[i],partNew1);
                    push!(partDetNew[i],partDetTemp[i][pci]);
                    push!(partNew[i],partNew2);
                    push!(partDetNew[i],partDetTemp[i][pci]);
                end
            end
        else
            for pci in 1:length(partCurrentTemp[i])
                push!(partNew[i],partCurrentTemp[i][pci]);
                push!(partDetNew[i],partDetTemp[i][pci]);
            end
        end
    end
    return partNew,partDetNew;
end

function combinePart(pData,disData,Ω,partCurrentTemp,partDetTemp,partHist)
    # aggregate some partition and anchor their values
end
