include("pullTightenBranch.jl");

# build the LP to obtain the longest path to each activity without crashing
function lpSolve(pData,iTarget)
    mp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0, CPX_PARAM_SIMDISPLAY = 0,CPX_PARAM_MIPDISPLAY = 0));
    @variable(mp, t[i in pData.II] >= 0);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]);
    @objective(mp, Min, t[iTarget]);
    solve(mp);

    return getobjectivevalue(mp);
end

# build the LP to obtain the longest path to each activity without crashing after the disruption
function lpSolveO1(pData,dDω,iTarget,brω)
    mp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0, CPX_PARAM_SIMDISPLAY = 0,CPX_PARAM_MIPDISPLAY = 0));
    @variable(mp, t[i in pData.II] >= 0);
    @constraint(mp, durationConstr1[k in pData.K; brω[findin(pData.II,k[1])[1]] >= 0], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]);
    @constraint(mp, tGbound[i in pData.II; brω[findin(pData.II,i)[1]] >= 0], t[i] >= dDω.H);
    @constraint(mp, durationConstr2[k in pData.K; brω[findin(pData.II,k[1])[1]] == -1], t[k[2]] - t[k[1]] >= pData.D[k[1]]);
    @objective(mp, Min, t[iTarget]);
    solve(mp);

    return getobjectivevalue(mp) - dDω.H;
end

function lpSolveO2(pData,dDω,iTarget,brω)
    mp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0, CPX_PARAM_SIMDISPLAY = 0,CPX_PARAM_MIPDISPLAY = 0));
    @variable(mp, t[i in pData.II] >= 0);
    @constraint(mp, tGbound[i in pData.II; brω[findin(pData.II,i)[1]] >= 0], t[i] >= dDω.H);
    @constraint(mp, durationConstr2[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]);
    @objective(mp, Min, t[iTarget]);
    solve(mp);

    return getobjectivevalue(mp) - dDω.H;
end

# bound tightening procedure
function buildTighten(pData,disData,Ω,cutSet,ub,lb,tmax,tmin,brInfo)
    bp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0, CPX_PARAM_SIMDISPLAY = 0,CPX_PARAM_MIPDISPLAY = 0));
    @variable(bp, tmin[i] <= t[i in pData.II] <= tmax[i]);
    @variable(bp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(bp, 0 <= F[i in pData.II, ω in Ω] <= 1);
    @variable(bp, 0 <= G[i in pData.II, ω in Ω] <= 1);
    @variable(bp, θ[ω in Ω] >= 0);

    @constraint(bp, upperbound, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω) <= ub);
    @constraint(bp, lowerbound, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω) >= lb);
    @constraint(bp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x[k[1],j]
        for j in pData.Ji[k[1]])));
    @constraint(bp, GFixed[i in pData.II,ω in Ω; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1],G[i,ω] == 1);
    @constraint(bp, FFixed[i in pData.II,ω in Ω; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1],F[i,ω] == 1);
    @constraint(bp, GPre[k in pData.K,ω in Ω], G[k[2],ω] >= G[k[1],ω]);
    @constraint(bp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(bp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(bp, tGnAnt[i in pData.II, ω in Ω; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], t[i] <= disData[ω].H + G[i,ω]*min(tmaxD[i] - disData[ω].H,M[i]));
    @constraint(bp, tFnAnt[i in pData.II, ω in Ω; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], t[i] >= disData[ω].H - F[i,ω]*(disData[ω].H - tminD[i]));
    @constraint(bp, FGnAnt[i in pData.II, ω in Ω; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], F[i,ω] + G[i,ω] == 1);
    # (πr,λr,γr,that,xhat,Ghat,vω)
    @constraint(bp, cuts[l in 1:length(cutSet),ω in Ω], θ[ω] >= cutSet[l][7][ω] + sum(cutSet[l][1][ω][i]*(t[i] - cutSet[l][4][i]) + cutSet[l][3][ω][i]*(G[i,ω] - cutSet[l][6][i,ω]) +
        sum(cutSet[l][2][ω][i,j]*(x[i,j] - cutSet[l][5][i,j]) for j in pData.Ji[i]) for i in pData.II));
    return bp;
end

function cutProc_Pull(pData,disData,Ω,TOL = 1e-6)
# initialize the first node
    nodeList = [];
    # add feasibility constraint up front
    M = Dict();
    for i in pData.II
        M[i] = lpSolve(pData,i);
    end
    brInfo = precludeRel(pData,disData,Ω);
    LB = 1e-6;
    UB = 999999;

    TOL = 1e-6;
    mIter = 500;

    tmaxD = Dict();
    tminD = Dict();
    maxH = maximum([disData[ω].H for ω in Ω]);
    for i in pData.II
        tmaxD[i] = maxH + M[i];
        tminD[i] = 0;
    end

    cutSet = [];
    nodeIni = nodeTypeP(LB,brInfo,tmaxD,tminD,cutSet,true);
    push!(nodeList,nodeIni);
    termCond = true;
    noI = 0;
    nCurrent = nodeList[1];

    bestUB = 999999;
    tbest = Dict();
    xbest = Dict();
    ωSeq,ωDict = getOmegaSeq(disData);

# start the iteration
    while ((bestUB - nCurrent.lbCost)/(nCurrent.lbCost) > TOL)&(termCond)
        noI += 1;
        # calculate the big M's for master/sub construction
        Mω1 = Dict();
        Mω2 = Dict();
        for ω in Ω
            # Mω[ω] = 10*M;
            dDω = disData[ω];
            brω = nCurrent.brInfo[:,ω];
            for i in pData.II
                Mω1[i,ω] = max(lpSolveO1(pData,dDω,i,brω),0);
                Mω2[i,ω] = max(lpSolveO2(pData,dDω,i,brω),0);
            end
        end

        # if the optimality gap is still large, continue the loop
        # construct the master/sub problem using the node information
        mp = Model(solver = CplexSolver(CPX_PARAM_EPAGAP = 1e-6,CPX_PARAM_EPRHS = 1e-9, CPX_PARAM_EPINT = 1e-9));
        @variable(mp, nCurrent.tminD[i] <= t[i in pData.II] <= nCurrent.tmaxD[i]);
        @variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
        @variable(mp, F[i in pData.II, ω in Ω], Bin);
        @variable(mp, G[i in pData.II, ω in Ω], Bin);
        @variable(mp, θ[ω in Ω] >= 0);

        @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

        @constraint(mp, objUB, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω) <= bestUB);
        @constraint(mp, objLB, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω) >= nCurrent.lbCost);

        @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x[k[1],j]
            for j in pData.Ji[k[1]])));
        @constraint(mp, GFixed1[i in pData.II,ω in Ω; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1],G[i,ω] == 1);
        @constraint(mp, GFixed2[i in pData.II,ω in Ω; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1],F[i,ω] == 1);
        @constraint(mp, GPre[k in pData.K,ω in Ω], G[k[2],ω] >= G[k[1],ω]);
        @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
        @constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
        @constraint(mp, tGnAnt[i in pData.II, ω in Ω; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], t[i] <= disData[ω].H + G[i,ω]*min(nCurrent.tmaxD[i] - disData[ω].H,M[i]));
        @constraint(mp, tFnAnt[i in pData.II, ω in Ω; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], t[i] >= disData[ω].H - F[i,ω]*(disData[ω].H - nCurrent.tminD[i]));
        @constraint(mp, FGnAnt[i in pData.II, ω in Ω], F[i,ω] + G[i,ω] == 1);
        @constraint(mp, cuts[l in 1:length(nCurrent.cutSet),ω in Ω], θ[ω] >= nCurrent.cutSet[l][7][ω] + sum(nCurrent.cutSet[l][1][ω][i]*(t[i] - nCurrent.cutSet[l][4][i]) + nCurrent.cutSet[l][3][ω][i]*(G[i,ω] - nCurrent.cutSet[l][6][i,ω]) +
            sum(nCurrent.cutSet[l][2][ω][i,j]*(x[i,j] - nCurrent.cutSet[l][5][i,j]) for j in pData.Ji[i]) for i in pData.II));

        xhatList = [];
        thatList = [];
        LBList = [];
        contBool = [true];
        cutNoList = [];

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

            UB = MathProgBase.cbgetobj(cb);
            noExplored = MathProgBase.cbgetexplorednodes(cb);
            push!(LBList,LB);
            println("--------------------------------------------------")
            println(statusCb," ",LB," ",UB);
            push!(xhatList,xhat);
            push!(thatList,that);

            # check the list contBool

            # if it converges to tolerance
            #if contBool[length(contBool)]
            if ((UB - LB)/LB > TOL)
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

                    @constraint(sp, tGbound1[i in pData.II; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1], tω[i] >= disData[ω].H);
                    @constraint(sp, tGbound2[i in pData.II; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1], tω[i] == getvalue(t[i]));
                    @constraint(sp, tGbound3[i in pData.II; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], tω[i] >= disData[ω].H*getvalue(G[i,ω]));
                    @constraint(sp, xGbound1[i in pData.II, j in pData.Ji[i]; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1], xω[i,j] == getvalue(x[i,j]));

                    @constraint(sp, xConstr[i in pData.II], sum(xω[i,j] for j in pData.Ji[i]) <= 1);
                    @constraint(sp, budgetConstr, sum(sum(xω[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);

                    @constraint(sp, tGnAnt1[i in pData.II; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], tω[i] <= getvalue(t[i]) + getvalue(G[i,ω])*Mω1[i,ω]);
                    @constraint(sp, tGnAnt2[i in pData.II; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], tω[i] >= getvalue(t[i]) - getvalue(G[i,ω])*Mω2[i,ω]);
                    @constraint(sp, xGnAnt1[i in pData.II, j in pData.Ji[i]; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], xω[i,j] <= getvalue(x[i,j]) + getvalue(G[i,ω]));
                    @constraint(sp, xGnAnt2[i in pData.II, j in pData.Ji[i]; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], xω[i,j] >= getvalue(x[i,j]) - getvalue(G[i,ω]));

                    @constraint(sp, sConstr1[i in pData.II, j in pData.Ji[i]; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], sω[i,j] - xω[i,j] <= 0);
                    @constraint(sp, sConstr2[i in pData.II, j in pData.Ji[i]; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], sω[i,j] <= getvalue(G[i,ω]));
                    @constraint(sp, sConstr3[i in pData.II, j in pData.Ji[i]; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0], sω[i,j] - xω[i,j] + 1 >= getvalue(G[i,ω]));
                    @constraint(sp, sConstr4[i in pData.II, j in pData.Ji[i]; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1], sω[i,j] == 0);
                    @constraint(sp, sConstr5[i in pData.II, j in pData.Ji[i]; nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1], sω[i,j] - xω[i,j] == 0);

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
                        if nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1
                            πr[ω][i] = getdual(tGbound2[i]);
                        elseif nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                            πr[ω][i] = getdual(tGnAnt1[i]) + getdual(tGnAnt2[i]);
                        else
                            πr[ω][i] = 0;
                        end
                        for j in pData.Ji[i]
                            # obtain λ's (for x)
                            if nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1
                                λr[ω][i,j] = getdual(xGbound1[i,j]);
                            elseif nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                                λr[ω][i,j] = getdual(xGnAnt1[i,j]) + getdual(xGnAnt2[i,j]);
                            else
                                λr[ω][i,j] = 0;
                            end
                        end
                        # obtain γ's (for G)
                        if i != 0
                            if nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1
                                γr[ω][i] = sum(disData[ω].d[i]*getdual(durationConstr)[k] for k in pData.K if k[1] == i);
                            elseif nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                                γr[ω][i] = sum(disData[ω].d[i]*getdual(durationConstr[k]) for k in pData.K if k[1] == i) +
                                    disData[ω].H*getdual(tGbound3[i]) + Mω1[i,ω]*(getdual(tGnAnt1)[i]) - Mω2[i,ω]*(getdual(tGnAnt2)[i]) +
                                    sum(getdual(xGnAnt1[i,j]) - getdual(xGnAnt2[i,j]) + getdual(sConstr2[i,j]) + getdual(sConstr3[i,j]) for j in pData.Ji[i]);
                            else
                                γr[ω][i] = 0;
                            end
                        else
                            if nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
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
                push!(nCurrent.cutSet,(πr,λr,γr,that,xhat,Ghat,vr));
            else
                return JuMP.StopTheSolver
            end

        end #innerCut

        function getInfo(cb)
            LB = MathProgBase.cbgetbestbound(cb);
            push!(LBList,LB);
            statusCb = MathProgBase.cbgetstate(cb);
            UB = MathProgBase.cbgetobj(cb);
            noExplored = MathProgBase.cbgetexplorednodes(cb);

            # after adding those cuts recalculate the brInfo and tbounds for every m iterations
            if (noExplored % mIter == 0)
                # push!(cutNoList,length(nCurrent.cutSet));
                # if length(cutNoList) == 1
                #     cutNoBool = true;
                # else
                #     cutNoBool = (cutNoList[length(cutNoList)] > cutNoList[length(cutNoList)-1]);
                # end
                # if cutNoBool
                # solve the new t bounds and brInfo
                println(UB," ",LB," ",length(nCurrent.cutSet)," ",noExplored);
                brInfoTemp = copy(nCurrent.brInfo);
                bp = buildTighten(pData,disData,Ω,nCurrent.cutSet,UB,LB,nCurrent.tmaxD,nCurrent.tminD,nCurrent.brInfo);
                maxt = Dict();
                mint = Dict();
                for i in pData.II
                    # change the objective function value
                    @objective(bp, Max, bp[:t][i]);
                    solve(bp);
                    #print(getvalue(bp[:t]));
                    maxt[i] = getobjectivevalue(bp);
                    @objective(bp, Min, bp[:t][i]);
                    solve(bp);
                    #print(getvalue(bp[:t]));
                    mint[i] = getobjectivevalue(bp);
                    for ω in Ω
                        if maxt[i] < disData[ω].H
                            if nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                                nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] = -1;
                            end
                        end
                        if mint[i] >= disData[ω].H
                            if nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                                nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] = 1;
                            end
                        end
                    end
                end
                nCurrent.brInfo = brInfoExt(pData,disData,Ω,nCurrent.brInfo,ωSeq);
                # if the brInfo and t bounds do not change then terminate
                contBoolT = false;
                for i in pData.II
                    #if the bounds are changed, then continue
                    for ω in Ω
                        if brInfoTemp[findin(pData.II,i)[1],findin(Ω,ω)[1]] != nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]]
                            contBoolT = true;
                        end
                    end
                    # println(i," ",nCurrent.tmaxD[i]," ",maxt[i]," ",nCurrent.tmaxD[i] > maxt[i]," ",nCurrent.tminD[i]," ",mint[i]," ",nCurrent.tminD[i] < mint[i]);
                    # if nCurrent.tmaxD[i] > maxt[i]
                    #     contBoolT = true;
                    #     nCurrent.tmaxD[i] = maxt[i];
                    # end
                    # if nCurrent.tminD[i] < mint[i]
                    #     contBoolT = true;
                    #     nCurrent.tminD[i] = mint[i];
                    # end
                end
                push!(contBool,contBoolT);
                println(contBool);
                lcB = length(contBool);
                if lcB > 3
                    contB = contBool[lcB];
                    for l in lcB:-1:lcB-3
                        contB = contB | contBool[l];
                    end
                else
                    contB = true;
                end
                if !contB
                    println("STOP!");
                    return JuMP.StopTheSolver
                end
            end
        end

        # solve the problem until the termination condition
        addlazycallback(mp,innerCut);
        addcutcallback(mp,getInfo);
        mpStatus = solve(mp);
        if mpStatus != :Infeasible
            # update the lower bound
            if LBList[length(LBList)] > nCurrent.lbCost
                nCurrent.lbCost = LBList[length(LBList)];
            end

            # if there is no solution it is okay and the process will still proceed
            tstar = getvalue(mp[:t]);
            xstar = getvalue(mp[:x]);
            UB = getobjectivevalue(mp);
            if bestUB > UB
                bestUB = UB;
                for i in pData.II
                    tbest[i] = tstar[i];
                    for j in pData.Ji[i]
                        xbest[i,j] = xstar[i,j];
                    end
                end
            end
            println("-------------------- Node $(noI) --------------------");
            println("Lower Bound: $(nCurrent.lbCost), Upper Bound: $(bestUB)");

            # branch to two nodes and update the nodeList
            if (mpStatus == :UserLimit)&(bestUB > nCurrent.lbCost)
                # branch if the problem is not solved to exact
                node1,node2 = branchPull(pData,disData,Ω,nCurrent,nCurrent.lbCost,bestUB);
                if node1.state
                    push!(nodeList,node1);
                end
                if node2.state
                    push!(nodeList,node2);
                end
            end
        end

        # select the next available node to probe
        shift!(nodeList);
        if nodeList == []
            termCond = false;
        else
            nCurrent = nodeList[1];
        end

    end

    # output the best solution and the best upper bound
    return tbest, xbest, bestUB;
end
