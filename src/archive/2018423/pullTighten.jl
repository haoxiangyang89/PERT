# pulling the binary variables to the first stage combined with bound tightening

function pullPartial(pData,disData,Ω,ϵ = 1e-6)
    # initialization
    M = Dict();
    for i in pData.II
        M[i] = lpSolve(pData,i);
    end
    # add feasibility constraint up front
    brInfo = precludeRel(pData,disData,Ω);

    LB = 1e-6;
    UB = 999999;
    TOL = 1e-6;
    tmaxD = Dict();
    tminD = Dict();
    maxH = maximum([disData[ω].H for ω in Ω]);
    for i in pData.II
        tmaxD[i] = maxH + M[i];
        tminD[i] = 0;
    end
    cutSetOri = [];

    while ((UB - LB)/LB > TOL)
        Mω1 = Dict();
        Mω2 = Dict();
        for ω in Ω
            # Mω[ω] = 10*M;
            dDω = disData[ω];
            brω = brInfo[:,ω];
            for i in pData.II
                Mω1[i,ω] = max(lpSolveO1(pData,dDω,i,brω),0);
                Mω2[i,ω] = max(lpSolveO2(pData,dDω,i,brω),0);
            end
        end
        # build the master problem
        #mp = Model(solver = GurobiSolver(Method = 0, OutputFlag = 0, IntFeasTol = 1e-9));
        mp = Model(solver = CplexSolver(CPX_PARAM_EPAGAP = 1e-6,CPX_PARAM_EPRHS = 1e-9, CPX_PARAM_EPINT = 1e-9));
        @variable(mp, tminD[i] <= t[i in pData.II] <= tmaxD[i]);
        @variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
        @variable(mp, F[i in pData.II, ω in Ω], Bin);
        @variable(mp, G[i in pData.II, ω in Ω], Bin);
        @variable(mp, θ[ω in Ω] >= 0);

        @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

        @constraint(mp, objUB, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω) <= UB);

        @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x[k[1],j]
            for j in pData.Ji[k[1]])));
        @constraint(mp, GFixed1[i in pData.II,ω in Ω; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1],G[i,ω] == 1);
        @constraint(mp, GFixed2[i in pData.II,ω in Ω; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1],F[i,ω] == 1);
        @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
        @constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    #    @constraint(mp, tGnAnt[i in pData.II, ω in Ω], t[i] <= (disData[ω].H - ϵ) + G[i,ω]*M);
        @constraint(mp, tGnAnt[i in pData.II, ω in Ω;brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], t[i] <= disData[ω].H + G[i,ω]*min(tmaxD[i] - disData[ω].H,M[i]));
        @constraint(mp, tFnAnt[i in pData.II, ω in Ω;brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], t[i] >= disData[ω].H - F[i,ω]*(disData[ω].H - tminD[i]));
        @constraint(mp, FGnAnt[i in pData.II, ω in Ω], F[i,ω] + G[i,ω] == 1);
        @constraint(mp, cuts[l in 1:length(cutSetOri),ω in Ω], θ[ω] >= cutSetOri[l][7][ω] + sum(cutSetOri[l][1][ω][i]*(t[i] - cutSetOri[l][4][i]) + cutSetOri[l][3][ω][i]*(G[i,ω] - cutSetOri[l][6][i,ω]) +
            sum(cutSetOri[l][2][ω][i,j]*(x[i,j] -cutSetOri[l][5][i,j]) for j in pData.Ji[i]) for i in pData.II));

        xhatList = [];
        thatList = [];
        GhatList = [];
        cutSet = [];
        UBList = [];
        LBList = [];
        nodesExploredList = [];

        function innerCut(cb)
            # callback function
            TOL = 1e-6;

            #if MathProgBase.cbgetbestbound(cb) > LB
            LB = MathProgBase.cbgetbestbound(cb);
            #end
            statusCb = MathProgBase.cbgetstate(cb);
            xhat = getvalue(x);
            that = getvalue(t);
            Ghat = getvalue(G);
            #UB = minimum(ubList);
            #if ubTemp < UB
            #    UB = ubTemp;
            #end
            UB = MathProgBase.cbgetobj(cb);
            noExplored = MathProgBase.cbgetexplorednodes(cb);
            push!(LBList,LB);
            push!(UBList,UB);
            println("--------------------------------------------------")
            println(statusCb," ",LB," ",UB);
            push!(xhatList,xhat);
            push!(thatList,that);
            push!(GhatList,Ghat);
            push!(nodesExploredList,noExplored);

            # if it converges to tolerance
            contBool = true;
            try
                impLB = (LBList[length(LBList)] - LBList[length(LBList) - 10])/(nodesExploredList[length(nodesExploredList)] - nodesExploredList[length(nodesExploredList) - 10]);
                if impLB <= 1e-8
                    contBool = false;
                end
            end
            if ((UB - LB)/LB > TOL)&(contBool)
                # set up the dual result dictionaries
                πr = Dict();
                λr = Dict();
                γr = Dict();
                vr = Dict();

                for ω in Ω
                    sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0, CPX_PARAM_SIMDISPLAY = 0,CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_LPMETHOD = 1));
                    @variable(sp, tω[i in pData.II] >= 0);
                    @variable(sp, 0 <= xω[i in pData.II, j in pData.Ji[i]] <= 1);
                    @variable(sp, 0 <= sω[i in pData.II, j in pData.Ji[i]] <= 1);

                    @constraint(sp, durationConstr[k in pData.K], tω[k[2]] - tω[k[1]] + sum(pData.D[k[1]]*pData.eff[k[1]][j]*xω[k[1],j] for j in pData.Ji[k[1]])
                        + sum(disData[ω].d[k[1]]*pData.eff[k[1]][j]*sω[k[1],j] for j in pData.Ji[k[1]])>= pData.D[k[1]] + disData[ω].d[k[1]]*getvalue(G[k[1],ω]));

                    @constraint(sp, tGbound1[i in pData.II; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1], tω[i] >= disData[ω].H);
                    @constraint(sp, tGbound2[i in pData.II; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1], tω[i] == getvalue(t[i]));
                    @constraint(sp, tGbound3[i in pData.II; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], tω[i] >= disData[ω].H*getvalue(G[i,ω]));
                    @constraint(sp, xGbound1[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1], xω[i,j] == getvalue(x[i,j]));

                    @constraint(sp, xConstr[i in pData.II], sum(xω[i,j] for j in pData.Ji[i]) <= 1);
                    @constraint(sp, budgetConstr, sum(sum(xω[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);

                    @constraint(sp, tGnAnt1[i in pData.II; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], tω[i] <= getvalue(t[i]) + getvalue(G[i,ω])*Mω1[i,ω]);
                    @constraint(sp, tGnAnt2[i in pData.II; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], tω[i] >= getvalue(t[i]) - getvalue(G[i,ω])*Mω2[i,ω]);
                    @constraint(sp, xGnAnt1[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], xω[i,j] <= getvalue(x[i,j]) + getvalue(G[i,ω]));
                    @constraint(sp, xGnAnt2[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], xω[i,j] >= getvalue(x[i,j]) - getvalue(G[i,ω]));

                    @constraint(sp, sConstr1[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], sω[i,j] - xω[i,j] <= 0);
                    @constraint(sp, sConstr2[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], sω[i,j] <= getvalue(G[i,ω]));
                    @constraint(sp, sConstr3[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], sω[i,j] - xω[i,j] + 1 >= getvalue(G[i,ω]));
                    @constraint(sp, sConstr4[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1], sω[i,j] == 0);
                    @constraint(sp, sConstr5[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1], sω[i,j] - xω[i,j] == 0);

                    @objective(sp, Min, tω[0]);
                    spStatus = solve(sp);
                    vω = getobjectivevalue(sp);
                    vr[ω] = vω;

                    # collect the dual variables
                    πr[ω] = Dict();
                    λr[ω] = Dict();
                    γr[ω] = Dict();
                    for i in pData.II
                        # obtain π's (for t)
                        if brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1
                            πr[ω][i] = getdual(tGbound2[i]);
                        elseif brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                            πr[ω][i] = getdual(tGnAnt1[i]) + getdual(tGnAnt2[i]);
                        else
                            πr[ω][i] = 0;
                        end
                        for j in pData.Ji[i]
                            # obtain λ's (for x)
                            if brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1
                                λr[ω][i,j] = getdual(xGbound1[i,j]);
                            elseif brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                                λr[ω][i,j] = getdual(xGnAnt1[i,j]) + getdual(xGnAnt2[i,j]);
                            else
                                λr[ω][i,j] = 0;
                            end
                        end
                        # obtain γ's (for G)
                        if i != 0
                            if brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1
                                γr[ω][i] = sum(disData[ω].d[i]*getdual(durationConstr)[k] for k in pData.K if k[1] == i);
                            elseif brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                                γr[ω][i] = sum(disData[ω].d[i]*getdual(durationConstr[k]) for k in pData.K if k[1] == i) +
                                    disData[ω].H*getdual(tGbound3[i]) + Mω1[i,ω]*(getdual(tGnAnt1)[i]) - Mω2[i,ω]*(getdual(tGnAnt2)[i]) +
                                    sum(getdual(xGnAnt1[i,j]) - getdual(xGnAnt2[i,j]) + getdual(sConstr2[i,j]) + getdual(sConstr3[i,j]) for j in pData.Ji[i]);
                            else
                                γr[ω][i] = 0;
                            end
                        else
                            if brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                                γr[ω][i] = disData[ω].H*getdual(tGbound3[i]) + Mω1[i,ω]*(getdual(tGnAnt1[i])) - Mω2[i,ω]*(getdual(tGnAnt2[i]));
                            else
                                γr[ω][i] = 0;
                            end
                        end
                    end
                    # @lazyconstraint(cb, θ[ω] >= vω + sum(πr[ω][i]*(t[i] - getvalue(t[i])) + γr[ω][i]*(G[i,ω] - getvalue(G[i,ω])) +
                    #     sum(λr[ω][i,j]*(x[i,j] - getvalue(x[i,j])) for j in pData.Ji[i]) for i in pData.II),localcut=true);
                    @lazyconstraint(cb, θ[ω] >= vω + sum(πr[ω][i]*(t[i] - getvalue(t[i])) + γr[ω][i]*(G[i,ω] - getvalue(G[i,ω])) +
                        sum(λr[ω][i,j]*(x[i,j] - getvalue(x[i,j])) for j in pData.Ji[i]) for i in pData.II));
                end
                push!(cutSet,(πr,λr,γr,that,xhat,Ghat,vr));
            else
                return JuMP.StopTheSolver;
            end

        end #innerCut

        # solve the master problem with the callback
        addlazycallback(mp,innerCut);
        statusMP = solve(mp);
        UB = min(UBList[length(UBList)],UB);
        LB = max(LBList[length(LBList)],LB);
        cutSetOri = append!(cutSetOri,cutSet);
        # set simple bounds on t and G
        if (UB - LB)/LB > TOL
            bp = buildTighten(pData,disData,Ω,cutSetOri,UB);
            maxt = Dict();
            mint = Dict();
            for i in pData.II
                # change the objective function value
                @objective(bp, Max, bp[:t][i]);
                solve(bp);
                maxt[i] = getobjectivevalue(bp);
                tmaxD[i] = maxt[i];
                @objective(bp, Min, bp[:t][i]);
                solve(bp);
                mint[i] = getobjectivevalue(bp);
                tminD[i] = mint[i];
                for ω in Ω
                    if maxt[i] < disData[ω].H
                        if brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                            brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] = -1;
                        end
                    end
                    if mint[i] >= disData[ω].H
                        if brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                            brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] = 1;
                        end
                    end
                end
            end
        end
    end # while

    bestObj = getobjectivevalue(mp);
    that = Dict();
    xhat = Dict();
    for i in pData.II
        that[i] = getvalue(mp[:t])[i];
        for j in pData.Ji[i]
            xhat[i,j] = getvalue(mp[:x])[i,j];
        end
    end

    return bestObj,that,xhat;

end
