# bound tightening procedure
function buildTighten(pData,disData,Ω,cutSet,ub,bigM = 999999)
    bp = Model(solver = GurobiSolver(OutputFlag = 0,Threads = 1));
    @variable(bp, t[i in pData.II] >= 0);
    @variable(bp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(bp, 0 <= F[i in pData.II, ω in Ω] <= 1);
    @variable(bp, 0 <= G[i in pData.II, ω in Ω] <= 1);
    @variable(bp, θ[ω in Ω] >= 0);

    @constraint(bp, upperbound, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω) <= ub);
    @constraint(bp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x[k[1],j]
        for j in pData.Ji[k[1]])));
    @constraint(bp, GFixed[i in pData.II,ω in Ω; brInfo[findfirst(x -> x == i, pData.II),findfirst(x -> x == ω, Ω)] == 1],G[i,ω] == 1);
    @constraint(bp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(bp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(bp, tGnAnt[i in pData.II, ω in Ω], t[i] <= disData[ω].H + G[i,ω]*bigM);
    @constraint(bp, tFnAnt[i in pData.II, ω in Ω], t[i] >= disData[ω].H - F[i,ω]*disData[ω].H);
    @constraint(bp, FGnAnt[i in pData.II, ω in Ω], F[i,ω] + G[i,ω] == 1);
    # (πr,λr,γr,that,xhat,Ghat,vω)
    @constraint(bp, cuts[l in 1:length(cutSet),ω in Ω], θ[ω] >= cutSet[l][7][ω] + sum(cutSet[l][1][ω][i]*(t[i] - cutSet[l][4][i]) + cutSet[l][3][ω][i]*(G[i,ω] - cutSet[l][6][i,ω]) +
        sum(cutSet[l][2][ω][i,j]*(x[i,j] -cutSet[l][5][i,j]) for j in pData.Ji[i]) for i in pData.II));
    return bp;
end

# pull the binary variables to the first stage
function pullDecomp(pData,disData,Ω,ϵ = 1e-6,bigM = 999999)
    # add feasibility constraint up front
    brInfo = precludeRel(pData,disData,Ω);

    # if the bigM is scenario generic, make it scenario specific
    if length(bigM) == 1
        bigMTemp = bigM;
        for i in pData.II
            for ω in Ω
                bigM[i,ω] = bigMTemp;
            end
        end
    end

    # build the master problem
    mp = Model(solver = GurobiSolver(Method = 0, IntFeasTol = 1e-9));
    #mp = Model(solver = CplexSolver(CPX_PARAM_EPAGAP = 1e-6,CPX_PARAM_EPRHS = 1e-9, CPX_PARAM_EPINT = 1e-9));
    @variable(mp, t[i in pData.II] >= 0);
    @variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(mp, G[i in pData.II, ω in Ω], Bin);
    @variable(mp, θ[ω in Ω] >= 0);

    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x[k[1],j]
        for j in pData.Ji[k[1]])));
    @constraint(mp, GFixed[i in pData.II,ω in Ω; brInfo[findfirst(x -> x == i, pData.II),findfirst(x -> x == ω, Ω)] == 1],G[i,ω] == 1);
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
#    @constraint(mp, tGnAnt[i in pData.II, ω in Ω], t[i] <= (disData[ω].H - ϵ) + G[i,ω]*M);
    @constraint(mp, tGnAnt[i in pData.II, ω in Ω], t[i] <= disData[ω].H + G[i,ω]*bigM);
    @constraint(mp, tFnAnt[i in pData.II, ω in Ω], t[i] >= disData[ω].H - (1-G[i,ω])*bigM);
    LB = -Inf;
    UB = Inf;

    xhatList = [];
    thatList = [];
    GhatList = [];
    ubList = [];
    cutSet = [];
    UBList = [];

    @everywhere function innerCut(cb)
        # callback function
        TOL = 1e-6;

        #if MathProgBase.cbgetbestbound(cb) > LB
        LB = MathProgBase.cbgetbestbound(cb);
        #end
        statusCb = MathProgBase.cbgetstate(cb);
        xhat = getvalue(x);
        that = getvalue(t);
        Ghat = getvalue(G);
        ubTemp = ubCalP(pData,disData,Ω,xhat,that,bigM);
        lbTemp = that[0]*pData.p0+sum(getvalue(θ[ω])*disData[ω].prDis for ω in Ω);
        push!(ubList,ubTemp);
        UB = minimum(ubList);
        if ubTemp < UB
           UB = ubTemp;
        end
        #UB = MathProgBase.cbgetobj(cb);
        push!(UBList,UB);
        println("--------------------------------------------------")
        println(statusCb," ",LB," ",UB," ",lbTemp," ",ubTemp);
        push!(xhatList,xhat);
        push!(thatList,that);
        push!(GhatList,Ghat);
        # if it converges to tolerance
        if (UB - LB)/LB > TOL
            # set up the dual result dictionaries
            πr = Dict();
            λr = Dict();
            γr = Dict();
            vr = Dict();

            for ω in Ω
                #sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0, CPX_PARAM_SIMDISPLAY = 0,CPX_PARAM_MIPDISPLAY = 0,CPX_PARAM_LPMETHOD = 1));
                sp = Model(solver = GurobiSolver(OutputFlag = 0));
                @variable(sp, tω[i in pData.II] >= 0);
                @variable(sp, 0 <= xω[i in pData.II, j in pData.Ji[i]] <= 1);
                @variable(sp, 0 <= sω[i in pData.II, j in pData.Ji[i]] <= 1);

                @constraint(sp, durationConstr[k in pData.K], tω[k[2]] - tω[k[1]] + sum(pData.D[k[1]]*pData.eff[k[1]][j]*xω[k[1],j] for j in pData.Ji[k[1]])
                    + sum(disData[ω].d[k[1]]*pData.eff[k[1]][j]*sω[k[1],j] for j in pData.Ji[k[1]])>= pData.D[k[1]] + disData[ω].d[k[1]]*getvalue(G[k[1],ω]));
                @constraint(sp, tGbound[i in pData.II], tω[i] >= disData[ω].H*Ghat[i,ω]);
                @constraint(sp, xConstr[i in pData.II], sum(xω[i,j] for j in pData.Ji[i]) <= 1);
                @constraint(sp, budgetConstr, sum(sum(xω[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
                @constraint(sp, tGnAnt1[i in pData.II], tω[i] <= that[i] + Ghat[i,ω]*bigM);
                @constraint(sp, tGnAnt2[i in pData.II], tω[i] >= that[i] - Ghat[i,ω]*bigM);
                @constraint(sp, xGnAnt1[i in pData.II, j in pData.Ji[i]], xω[i,j] <= xhat[i,j] + Ghat[i,ω]);
                @constraint(sp, xGnAnt2[i in pData.II, j in pData.Ji[i]], xω[i,j] >= xhat[i,j] - Ghat[i,ω]);
                @constraint(sp, sConstr1[i in pData.II, j in pData.Ji[i]], sω[i,j] - xω[i,j] <= 0);
                @constraint(sp, sConstr2[i in pData.II, j in pData.Ji[i]], sω[i,j] <= Ghat[i,ω]);
                @constraint(sp, sConstr3[i in pData.II, j in pData.Ji[i]], sω[i,j] - xω[i,j] + 1 >= Ghat[i,ω]);

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
                    πr[ω][i] = getdual(tGnAnt1)[i] + getdual(tGnAnt2)[i];
                    for j in pData.Ji[i]
                        # obtain λ's (for x)
                        λr[ω][i,j] = getdual(xGnAnt1)[i,j] + getdual(xGnAnt2)[i,j];
                    end
                    # obtain γ's (for G)
                    if i != 0
                        γr[ω][i] = sum(disData[ω].d[i]*getdual(durationConstr)[k] for k in pData.K if k[1] == i) +
                            disData[ω].H*getdual(tGbound)[i] + bigM*(getdual(tGnAnt1)[i]) - bigM*(getdual(tGnAnt2)[i]) +
                            sum(getdual(xGnAnt1)[i,j] - getdual(xGnAnt2)[i,j] + getdual(sConstr2)[i,j] + getdual(sConstr3)[i,j] for j in pData.Ji[i]);
                    else
                        γr[ω][i] = disData[ω].H*getdual(tGbound)[i] + bigM*(getdual(tGnAnt1)[i]) - bigM*(getdual(tGnAnt2)[i]);
                    end
                end
                # @lazyconstraint(cb, θ[ω] >= vω + sum(πr[ω][i]*(t[i] - that[i]) + γr[ω][i]*(G[i,ω] - Ghat[i,ω]) +
                #     sum(λr[ω][i,j]*(x[i,j] - xhat[i,j]) for j in pData.Ji[i]) for i in pData.II),localcut=true);
                @lazyconstraint(cb, θ[ω] >= vω + sum(πr[ω][i]*(t[i] - that[i]) + γr[ω][i]*(G[i,ω] - Ghat[i,ω]) +
                    sum(λr[ω][i,j]*(x[i,j] - xhat[i,j]) for j in pData.Ji[i]) for i in pData.II));
            end
            push!(cutSet,(πr,λr,γr,that,xhat,Ghat,vr));
            # set simple bounds on t and G
            bp = buildTighten(pData,disData,Ω,cutSet,UB);
            maxt = Dict();
            mint = Dict();
            for i in pData.II
                # change the objective function value
                @objective(bp, Max, bp[:t][i]);
                solve(bp);
                maxt[i] = getobjectivevalue(bp);
                @objective(bp, Min, bp[:t][i]);
                solve(bp);
                mint[i] = getobjectivevalue(bp);
                for ω in Ω
                    if maxt[i] < disData[ω].H
                        if brInfo[findfirst(x -> x == i, pData.II),findfirst(x -> x == ω, Ω)] == 0
                            brInfo[findfirst(x -> x == i, pData.II),findfirst(x -> x == ω, Ω)] = -1;
                            @lazyconstraint(cb, G[i,ω] == 0);
                        end
                    end
                    if mint[i] >= disData[ω].H
                        if brInfo[findfirst(x -> x == i, pData.II),findfirst(x -> x == ω, Ω)] == 0
                            brInfo[findfirst(x -> x == i, pData.II),findfirst(x -> x == ω, Ω)] = 1;
                            @lazyconstraint(cb, G[i,ω] == 1);
                        end
                    end
                end
            end
        else
            return JuMP.StopTheSolver;
        end
    end #innerCut

    # solve the master problem with the callback
    addlazycallback(mp,innerCut);
    solve(mp);

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
