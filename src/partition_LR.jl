function lapConstruct(pData,disData,Ω,cutSet,partCurrentTemp,partDet,μ,M = 999999,FPre = [],intInd = 0)
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

    lap = Model(solver = GurobiSolver(OutputFlag = 0));

    @variable(lap, t[i in pData.II] >= 0);
    @variable(lap, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
    if intInd == 0
        @variable(lap, 0 <= G[i in pData.II, ω in Ω] <= 1);
    else
        @variable(lap, G[i in pData.II, ω in Ω], Bin);
    end
    @variable(lap, F[i in pData.II, m in 1:partNo[i]],Bin);
    if FPre != []
        @constraint(lap, Fbind[i in pData.II,τ in 1:length(partCurrentTemp[i])], F[i,τ] == FPre[i,τ]);
    end

    @variable(lap, θ[ω in Ω] >= 0);

    @objective(lap, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] + sum(μ[i,ω]*(G[i,ω] - F[i,partRev[i][ω]]) for i in pData.II) for ω in Ω));

    @constraint(lap, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(lap, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(lap, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    # logic constraints between G and t
    @constraint(lap, tGcons1[i in pData.II, ω in Ω], t[i] <= disData[ω].H + G[i,ω]*M);
    @constraint(lap, tGcons2[i in pData.II, ω in Ω], t[i] >= disData[ω].H - (1 - G[i,ω])*M);

    # logic constraints of G
    @constraint(lap, GijCon[i in pData.II, j in pData.Succ[i],ω in Ω], G[i,ω] <= G[j,ω]);
    @constraint(lap, Gω12Con[i in pData.II, ω in 1:(length(Ω) - 1)], G[i,ω] >= G[i,ω + 1]);
    @constraint(lap, Gdet1[i in pData.II,ω in Ω; partDet[i][partRev[i][ω]] == 1], G[i,ω] == 1);
    #@constraint(lap, Fdet1[i in pData.II,m in 1:partNo[i]; partDet[i][m] == 1], F[i,m] == 1);
    #@constraint(lap, Fω12Con[i in pData.II, m in 1:(partNo[i]-1)], F[i,m] >= F[i,m + 1]);

    # add cuts
    @constraint(lap, cuts[ω in Ω, nc in 1:length(cutSet[ω])], θ[ω] >= cutSet[ω][nc][4] +
        sum(cutSet[ω][nc][1][i]*(lap[:t][i] - cutSet[ω][nc][5][i]) for i in pData.II) +
        sum(sum(cutSet[ω][nc][3][i,j]*(lap[:x][i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II) +
        sum(cutSet[ω][nc][2][i]*(lap[:G][i,ω] - cutSet[ω][nc][7][i]) for i in pData.II));

    return lap;
end

function lapConstruct2(pData,disData,Ω,cutSet,partCurrentTemp,partDet,M = 999999,FPre = [])
    # dual model to solve the Lagrangian problem analytically
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
    lap = Model(solver = GurobiSolver());

    @variable(lap, λD[k in pData.K] >= 0);
    @variable(lap, λB <= 0);
    @variable(lap, λxbar[i in pData.II] <= 0);
    @variable(lap, λtG1[i in pData.II, ω in Ω] <= 0);
    @variable(lap, λtG2[i in pData.II, ω in Ω] >= 0);
    @variable(lap, λc[ω in Ω, n in 1:length(cutSet[ω])] >= 0);
    @variable(lap, λGbar[i in pData.II, ω in Ω] <= 0);
    @variable(lap, λG1[i in pData.II, ω in Ω] >= 0);
    @variable(lap, λG2[i in pData.II, j in pData.Succ[i], ω in Ω] <= 0);
    @variable(lap, λG3[i in pData.II, ω in Ω]);
    @variable(lap, μ[i in pData.II, ω in Ω]);
    @variable(lap, δ[i in pData.II, τ in 1:length(partCurrentTemp[i])] <= 0);

    # G problem
    @constraint(lap, t0Constr, pData.p0 - sum(λD[(j,0)] for j in pData.Pre[0]) - sum(λtG1[0,ω] + λtG2[0,ω] for ω in Ω) +
        sum(sum(λc[ω,n]*cutSet[ω][n][1][0] for n in 1:length(cutSet[ω])) for ω in Ω) >= 0);
    @constraint(lap, tiConstr[i in pData.II; i != 0], sum(λD[(i,j)] for j in pData.Succ[i]) - sum(λD[(j,i)] for j in pData.Pre[i]) -
        sum(λtG1[i,ω] + λtG2[i,ω] for ω in Ω) + sum(sum(λc[ω,n]*cutSet[ω][n][1][i] for n in 1:length(cutSet[ω])) for ω in Ω) >= 0);
    @constraint(lap, xConstr[i in pData.II, j in pData.Ji[i]], -λxbar[i] - λB - sum(λD[(i,k)]*pData.D[i]*pData.eff[i][j] for k in pData.Succ[i]) +
        sum(sum(λc[ω,n]*cutSet[ω][n][3][i,j] for n in 1:length(cutSet[ω])) for ω in Ω) >= 0);
    @constraint(lap, θConstr[ω in Ω], disData[ω].prDis - sum(λc[ω,n] for n in 1:length(cutSet[ω])) >= 0);

    coeffPart1 = Dict();
    for i in pData.II
        for ω in Ω
            if partDet[i][partRev[i][ω]] == 1
                coeffPart1[i,ω] = -1;
            else
                coeffPart1[i,ω] = 0;
            end
        end
    end
    @constraint(lap, GConstr[i in pData.II, ω in Ω;(ω != 1)&(ω != Ω[length(Ω)])], μ[i,ω] + M*(λtG1[i,ω] + λtG2[i,ω]) + sum(cutSet[ω][n][2][i]*λc[ω,n] for n in 1:length(cutSet[ω])) -
        λGbar[i,ω] + λG1[i,ω] - λG1[i,ω + 1] - sum(λG2[i,j,ω] for j in pData.Succ[i]) + sum(λG2[j,i,ω] for j in pData.Pre[i]) + coeffPart1[i,ω]*λG3[i,ω] >= 0);
    @constraint(lap, GConstr1[i in pData.II], μ[i,1] + M*(λtG1[i,1] + λtG2[i,1]) + sum(cutSet[1][n][2][i]*λc[1,n] for n in 1:length(cutSet[1])) -
        λGbar[i,1] - λG1[i,2] - sum(λG2[i,j,1] for j in pData.Succ[i]) + sum(λG2[j,i,1] for j in pData.Pre[i]) + coeffPart1[i,1]*λG3[i,1] >= 0);
    @constraint(lap, GConstrΩ[i in pData.II], μ[i,Ω[length(Ω)]] + M*(λtG1[i,Ω[length(Ω)]] + λtG2[i,Ω[length(Ω)]]) + sum(cutSet[Ω[length(Ω)]][n][2][i]*λc[Ω[length(Ω)],n] for n in 1:length(cutSet[Ω[length(Ω)]])) -
        λGbar[i,Ω[length(Ω)]] + λG1[i,Ω[length(Ω)]] - sum(λG2[i,j,Ω[length(Ω)]] for j in pData.Succ[i]) + sum(λG2[j,i,Ω[length(Ω)]] for j in pData.Pre[i]) + coeffPart1[i,Ω[length(Ω)]]*λG3[i,Ω[length(Ω)]] >= 0);

    # F problem
    @constraint(lap, δμConstr[i in pData.II, τ in 1:length(partCurrentTemp[i])], δ[i,τ] <= -sum(μ[i,ω] for ω in partCurrentTemp[i][τ]));
    if FPre == []
        @constraint(lap, δConstr1[i in pData.II, τ in 1:length(partCurrentTemp[i]);partDet[i][τ] == 1], δ[i,τ] == -sum(μ[i,ω] for ω in partCurrentTemp[i][τ]));
    else
        @constraint(lap, δConstr1[i in pData.II, τ in 1:length(partCurrentTemp[i])], δ[i,τ] == -FPre[i,τ]*sum(μ[i,ω] for ω in partCurrentTemp[i][τ]));
    end

    @objective(lap, Max, sum(λD[k]*pData.D[k[1]] for k in pData.K) + pData.B*λB + sum(λxbar[i] for i in pData.II) +
        sum(sum(disData[ω].H*λtG1[i,ω] + (disData[ω].H - M)*λtG2[i,ω] + λGbar[i,ω] for ω in Ω) for i in pData.II) +
        sum(sum((cutSet[ω][n][4] - sum(cutSet[ω][n][1][i]*cutSet[ω][n][5][i] for i in pData.II) -
            sum(sum(cutSet[ω][n][3][i,j]*cutSet[ω][n][6][i,j] for j in pData.Ji[i]) for i in pData.II) -
            sum(cutSet[ω][n][2][i]*cutSet[ω][n][7][i] for i in pData.II))*λc[ω,n] for n in 1:length(cutSet[ω])) for ω in Ω) +
        sum(sum(δ[i,τ] for τ in 1:length(partCurrentTemp[i])) for i in pData.II) +
        sum(sum(sum(λG3[i,ω] for ω in partCurrentTemp[i][τ]) for τ in 1:length(partCurrentTemp[i]) if partDet[i][τ] == 1) for i in pData.II));

    solve(lap);
    μp = Dict();
    for ω in Ω
        for i in pData.II
            μp[i,ω] = getvalue(lap[:μ][i,ω]);
        end
    end
    return μp;
end

function lapConstructFix(pData,disData,Ω,cutSet,partCurrentTemp,partDet,M = 999999,FPre = [])
    # dual model to solve the Lagrangian problem analytically
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
    lap = Model(solver = GurobiSolver());

    @variable(lap, λD[k in pData.K] >= 0);
    @variable(lap, λB <= 0);
    @variable(lap, λxbar[i in pData.II] <= 0);
    @variable(lap, λtG1[i in pData.II, ω in Ω] <= 0);
    @variable(lap, λtG2[i in pData.II, ω in Ω] >= 0);
    @variable(lap, λc[ω in Ω, n in 1:length(cutSet[ω])] >= 0);
    @variable(lap, λGbar[i in pData.II, ω in Ω] <= 0);
    @variable(lap, λG1[i in pData.II, ω in Ω] >= 0);
    @variable(lap, λG2[i in pData.II, j in pData.Succ[i], ω in Ω] <= 0);
    @variable(lap, λG3[i in pData.II, ω in Ω]);
    @variable(lap, μ[i in pData.II, ω in Ω]);
    @variable(lap, δ[i in pData.II, τ in 1:length(partCurrentTemp[i])] <= 0);

    # G problem
    @constraint(lap, t0Constr, pData.p0 - sum(λD[(j,0)] for j in pData.Pre[0]) - sum(λtG1[0,ω] + λtG2[0,ω] for ω in Ω) +
        sum(sum(λc[ω,n]*cutSet[ω][n][1][0] for n in 1:length(cutSet[ω])) for ω in Ω) >= 0);
    @constraint(lap, tiConstr[i in pData.II; i != 0], sum(λD[(i,j)] for j in pData.Succ[i]) - sum(λD[(j,i)] for j in pData.Pre[i]) -
        sum(λtG1[i,ω] + λtG2[i,ω] for ω in Ω) + sum(sum(λc[ω,n]*cutSet[ω][n][1][i] for n in 1:length(cutSet[ω])) for ω in Ω) >= 0);
    @constraint(lap, xConstr[i in pData.II, j in pData.Ji[i]], -λxbar[i] - λB - sum(λD[(i,k)]*pData.D[i]*pData.eff[i][j] for k in pData.Succ[i]) +
        sum(sum(λc[ω,n]*cutSet[ω][n][3][i,j] for n in 1:length(cutSet[ω])) for ω in Ω) >= 0);
    @constraint(lap, θConstr[ω in Ω], disData[ω].prDis - sum(λc[ω,n] for n in 1:length(cutSet[ω])) >= 0);

    coeffPart1 = Dict();
    for i in pData.II
        for ω in Ω
            if partDet[i][partRev[i][ω]] == 1
                coeffPart1[i,ω] = -1;
            else
                coeffPart1[i,ω] = 0;
            end
        end
    end
    @constraint(lap, GConstr[i in pData.II, ω in Ω;(ω != 1)&(ω != Ω[length(Ω)])], μ[i,ω] + M*(λtG1[i,ω] + λtG2[i,ω]) + sum(cutSet[ω][n][2][i]*λc[ω,n] for n in 1:length(cutSet[ω])) -
        λGbar[i,ω] + λG1[i,ω] - λG1[i,ω + 1] - sum(λG2[i,j,ω] for j in pData.Succ[i]) + sum(λG2[j,i,ω] for j in pData.Pre[i]) + coeffPart1[i,ω]*λG3[i,ω] >= 0);
    @constraint(lap, GConstr1[i in pData.II], μ[i,1] + M*(λtG1[i,1] + λtG2[i,1]) + sum(cutSet[1][n][2][i]*λc[1,n] for n in 1:length(cutSet[1])) -
        λGbar[i,1] - λG1[i,2] - sum(λG2[i,j,1] for j in pData.Succ[i]) + sum(λG2[j,i,1] for j in pData.Pre[i]) + coeffPart1[i,1]*λG3[i,1] >= 0);
    @constraint(lap, GConstrΩ[i in pData.II], μ[i,Ω[length(Ω)]] + M*(λtG1[i,Ω[length(Ω)]] + λtG2[i,Ω[length(Ω)]]) + sum(cutSet[Ω[length(Ω)]][n][2][i]*λc[Ω[length(Ω)],n] for n in 1:length(cutSet[Ω[length(Ω)]])) -
        λGbar[i,Ω[length(Ω)]] + λG1[i,Ω[length(Ω)]] - sum(λG2[i,j,Ω[length(Ω)]] for j in pData.Succ[i]) + sum(λG2[j,i,Ω[length(Ω)]] for j in pData.Pre[i]) + coeffPart1[i,Ω[length(Ω)]]*λG3[i,Ω[length(Ω)]] >= 0);

    # added μ monotonicity
    @constraint(lap, μConstr1[i in pData.II, τ in 1:length(partCurrentTemp[i]),ω in partCurrentTemp[i][τ];ω != partCurrentTemp[i][τ][length(partCurrentTemp[i][τ])]], μ[i,ω] >= μ[i,ω+1]);
    #@constraint(lap, μConstr2[i in pData.II, j in pData.Succ[i], ω in Ω], μ[i,ω] >= μ[j,ω]);

    # F problem
    @constraint(lap, δμConstr[i in pData.II, τ in 1:length(partCurrentTemp[i])], δ[i,τ] <= -sum(μ[i,ω] for ω in partCurrentTemp[i][τ]));
    if FPre == []
        @constraint(lap, δConstr1[i in pData.II, τ in 1:length(partCurrentTemp[i]);partDet[i][τ] == 1], δ[i,τ] == -sum(μ[i,ω] for ω in partCurrentTemp[i][τ]));
    else
        @constraint(lap, δConstr1[i in pData.II, τ in 1:length(partCurrentTemp[i])], δ[i,τ] == -FPre[i,τ]*sum(μ[i,ω] for ω in partCurrentTemp[i][τ]));
    end

    @objective(lap, Max, sum(λD[k]*pData.D[k[1]] for k in pData.K) + pData.B*λB + sum(λxbar[i] for i in pData.II) +
        sum(sum(disData[ω].H*λtG1[i,ω] + (disData[ω].H - M)*λtG2[i,ω] + λGbar[i,ω] for ω in Ω) for i in pData.II) +
        sum(sum((cutSet[ω][n][4] - sum(cutSet[ω][n][1][i]*cutSet[ω][n][5][i] for i in pData.II) -
            sum(sum(cutSet[ω][n][3][i,j]*cutSet[ω][n][6][i,j] for j in pData.Ji[i]) for i in pData.II) -
            sum(cutSet[ω][n][2][i]*cutSet[ω][n][7][i] for i in pData.II))*λc[ω,n] for n in 1:length(cutSet[ω])) for ω in Ω) +
        sum(sum(δ[i,τ] for τ in 1:length(partCurrentTemp[i])) for i in pData.II) +
        sum(sum(sum(λG3[i,ω] for ω in partCurrentTemp[i][τ]) for τ in 1:length(partCurrentTemp[i]) if partDet[i][τ] == 1) for i in pData.II));

    solve(lap);
    μp = Dict();
    for ω in Ω
        for i in pData.II
            μp[i,ω] = getvalue(lap[:μ][i,ω]);
        end
    end
    return μp;
end

function lapSub(pData,disData,Ω,cutSet,partCurrentTemp,partDet,μ,M = 999999)
    # dual model to solve the Lagrangian problem analytically
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

    sLap = Model(solver = GurobiSolver());
    @variables(sLap, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
      0 <= G[i in pData.II,ω in Ω] <= 1
    end);
    @constraint(sLap, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sLap, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(sLap, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    # logic constraints between G and t
    @constraint(sLap, tGcons1[i in pData.II, ω in Ω], t[i] <= disData[ω].H + G[i,ω]*M);
    @constraint(sLap, tGcons2[i in pData.II, ω in Ω], t[i] >= disData[ω].H - (1 - G[i,ω])*M);

    # logic constraints of G
    @constraint(sLap, GijCon[i in pData.II, j in pData.Succ[i],ω in Ω], G[i,ω] <= G[j,ω]);
    @constraint(sLap, Gω12Con[i in pData.II, ω in 1:(length(Ω) - 1)], G[i,ω] >= G[i,ω + 1]);
    @constraint(sLap, Gdet1[i in pData.II,ω in Ω; partDet[i][partRev[i][ω]] == 1], G[i,ω] == 1);
    #@constraint(sLap, Gdet0[i in pData.II,m in 1:partNo[i]; partDet[i][m] == -1], G[i,m] == 0);

    @constraint(sLap, cuts[ω in Ω, nc in 1:length(cutSet[ω])], θ[ω] >= cutSet[ω][nc][4] +
        sum(cutSet[ω][nc][1][i]*(t[i] - cutSet[ω][nc][5][i]) for i in pData.II) +
        sum(sum(cutSet[ω][nc][3][i,j]*(x[i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II) +
        sum(cutSet[ω][nc][2][i]*(G[i,ω] - cutSet[ω][nc][7][i]) for i in pData.II));

    @objective(sLap, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] + sum(μ[i,ω]*G[i,ω] for i in pData.II) for ω in Ω));
    return sLap;
end

function lapConstructCut(pData,disData,Ω,cutSet,partCurrentTemp,partDet,ub,M = 999999,FPre = [])
    mLap = Model(solver = GurobiSolver());
    @variable(mLap, -M <= μ[i in pData.II, ω in Ω] <= M);
    @variable(mLap, fF[i in pData.II, τ in 1:length(partCurrentTemp[i])]);
    @variable(mLap, fG <= 2*ub);

    if FPre == []
        @constraint(mLap, fF1[i in pData.II, τ in 1:length(partCurrentTemp[i])], fF[i,τ] <= 0);
        @constraint(mLap, fF2[i in pData.II, τ in 1:length(partCurrentTemp[i])], fF[i,τ] <= -sum(μ[i,ω] for ω in partCurrentTemp[i][τ]));
    else
        @constraint(mLap, fF2[i in pData.II, τ in 1:length(partCurrentTemp[i])], fF[i,τ] == -FPre[i,τ]*sum(μ[i,ω] for ω in partCurrentTemp[i][τ]));
    end
    @constraint(mLap, μConstr1[i in pData.II, τ in 1:length(partCurrentTemp[i]),ω in partCurrentTemp[i][τ];ω != partCurrentTemp[i][τ][length(partCurrentTemp[i][τ])]], μ[i,ω] >= μ[i,ω+1]);
    #@constraint(mLap, μConstr2[i in pData.II, j in pData.Succ[i], ω in Ω], μ[i,ω] >= μ[j,ω]);

    @objective(mLap, Max, fG + sum(sum(fF[i,τ] for τ in 1:length(partCurrentTemp[i])) for i in pData.II));

    stopBool = true;
    while stopBool
        solve(mLap);
        ubM = getobjectivevalue(mLap);
        μhat = Dict();
        for i in pData.II
            for ω in Ω
                μhat[i,ω] = getvalue(mLap[:μ][i,ω]);
            end
        end

        # generate the cut for fG
        sLap = lapSub(pData,disData,Ω,cutSet,partCurrentTemp,partDet,μhat,400);
        solve(sLap);
        sLapV = getobjectivevalue(sLap);
        if FPre == []
            lbIter = sLapV + sum(sum(min(0,sum(-μhat[i,ω] for ω in partCurrentTemp[i][τ])) for τ in 1:length(partCurrentTemp[i])) for i in pData.II);
        else
            lbIter = sLapV + sum(sum(-FPre[i,τ]*sum(μhat[i,ω] for ω in partCurrentTemp[i][τ]) for τ in 1:length(partCurrentTemp[i])) for i in pData.II);
        end
        Ghat = Dict();
        for i in pData.II
            for ω in Ω
                Ghat[i,ω] = getvalue(sLap[:G][i,ω]);
            end
        end

        # append the cut to mLap
        if ubM - lbIter < 1e-5
            stopBool = false;
        else
            @constraint(mLap, fG <= sLapV + sum(sum(Ghat[i,ω]*(mLap[:μ][i,ω] - μhat[i,ω]) for ω in Ω) for i in pData.II));
        end
        # @objective(mLap, Max, fG + sum(sum(fF[i,τ] for τ in 1:length(partCurrentTemp[i])) for i in pData.II) -
        #     1*(sum(sum((μ[i,ω] - μhat[i,ω])^2 for ω in Ω) for i in pData.II)));
    end
    μhat = Dict();
    for i in pData.II
        for ω in Ω
            μhat[i,ω] = getvalue(mLap[:μ][i,ω]);
        end
    end
    return μhat,mLap;
end

function lagrangianPart(pData,disData,Ω,cutSet,partCurrentTemp,partDet,ub,M = 999999)
    # solving the Lagrangian relaxation to obtain the dual evaluation
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

    # initialize a μ0
    μ = Dict();
    for ω in Ω
        for i in pData.II
            μ[i,ω] = 0;
        end
    end
    stopBool = false;
    β = 0.1;
    μList = [μ];
    lbList = [0.0];
    lbCounter = 0;
    iterCounter = 0;
    μbest = [];
    while !stopBool
        # construct Lagrangian relaxation and solve for the subgradient of μ
        lap = lapConstruct(pData,disData,Ω,cutSet,partCurrentTemp,partDet,μ,400);
        lapStatus = solve(lap);
        lb = getobjectivevalue(lap);
        Ginc = getvalue(lap[:G]);
        Finc = getvalue(lap[:F]);
        # obtain the inner problem solution and update μ
        deno = sum(sum((Ginc[i,ω] - Finc[i,partRev[i][ω]])^2 for i in pData.II) for ω in Ω);
        if lb > maximum(lbList)
            push!(lbList,lb);
            lbCounter = 0;
            μbest = copy(μ);
        else
            push!(lbList,lb);
            lbCounter += 1;
        end
        if lbCounter >= 20
            β = β/2;
            lbCounter = 0;
        end
        # update the dual price
        α = β*(ub - lb)/(deno);
        for ω in Ω
            for i in pData.II
                μ[i,ω] = μ[i,ω] + α*(Ginc[i,ω] - Finc[i,partRev[i][ω]]);
            end
        end
        # check if the stopping criterion is met
        iterCounter += 1;
        if iterCounter > 120
            stopBool = true;
        end
    end

    # split the partition
    partNew = Dict();
    partDetNew = Dict();
    for i in pData.II
        partNew[i] = [];
        partDetNew[i] = [];
        for pci in 1:length(partCurrentTemp[i])
            push!(partNew[i],partCurrentTemp[i][pci]);
            push!(partDetNew[i],partDet[i][pci]);
        end
    end
    for i in pData.II
        for m in 1:partNo[i]
            # find the partition to break
            if partDet[i][m] == 0
                # find the location to break, should always go from negative to positive
                stopPt = 0;
                for ω in partCurrentTemp[i][m]
                    if μ[i,ω] > 0
                        stopPt += 1;
                    end
                end
                if (stopPt != 0)|(stopPt != length(partCurrentTemp[i][m]))
                    # there is a beneficial spliting point
                    splice!(partNew[i],m,[partNew[i][m][1:stopPt],partNew[i][m][(stopPt + 1):length(partNew[i][m])]]);
                    splice!(partDetNew[i],m,[0,0]);
                end
            end
        end
    end
end
