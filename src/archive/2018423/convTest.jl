# test the convergence of the Benders decomposition
function subIntG(pData,dDω,xhat,that,Ghatω)
    # solve the MIP recourse problem
    M = sum(max(pData.D[i],pData.D[i]+dDω.d[i]) for i in pData.II if i != 0);

    # sp = Model(solver = GurobiSolver(OutputFlag = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);

    # add the basic sub problem constraints
    @constraint(sp, fixT[i in pData.II; Ghatω[i] == 0], t[i] == that[i]);
    @constraint(sp, boundT[i in pData.II; Ghatω[i] == 1], t[i] >= dDω.H);
    @constraint(sp, fixX[i in pData.II,j in pData.Ji[i]; Ghatω[i] == 0], x[i,j] == xhat[i,j]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K; Ghatω[k[1]] == 0], t[k[2]] - t[k[1]] >= pData.D[k[1]]*
        (1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(sp, durationConstr2[k in pData.K; Ghatω[k[1]] == 1], t[k[2]] - t[k[1]] >= (pData.D[k[1]]+dDω.d[k[1]])*
        (1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));

    @objective(sp, Min, t[0]);
    solve(sp);
    return getobjectivevalue(sp);
end

function ubCalG(pData,disData,Ω,xhat,that,Ghat)
    # calculate the upper bound of the problem given the master solution
    ubCost = that[0]*pData.p0;
    for ω in Ω
        cω = subIntG(pData,disData[ω],xhat,that,Ghat[ω]);
        ubCost += disData[ω].prDis*cω;
    end
    return ubCost;
end

M = Dict();
for i in pData.II
    M[i] = lpSolve(pData,i);
end
# add feasibility constraint up front
brInfo = precludeRel(pData,disData,Ω);

# build the master problem
#mp = Model(solver = GurobiSolver(Method = 0, OutputFlag = 0, IntFeasTol = 1e-9));
mp = Model(solver = CplexSolver(CPX_PARAM_EPAGAP = 1e-6,CPX_PARAM_EPRHS = 1e-9, CPX_PARAM_EPINT = 1e-9,CPX_PARAM_SCRIND = 0, CPX_PARAM_SIMDISPLAY = 0,CPX_PARAM_MIPDISPLAY = 0));
@variable(mp, t[i in pData.II] >= 0);
@variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
@variable(mp, F[i in pData.II, ω in Ω], Bin);
@variable(mp, G[i in pData.II, ω in Ω], Bin);
@variable(mp, θ[ω in Ω] >= 0);
# @variable(mp, θ >= 0);

@objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));
# @objective(mp, Min, pData.p0*t[0] + θ);

@constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x[k[1],j]
    for j in pData.Ji[k[1]])));
@constraint(mp, GFixed[i in pData.II,ω in Ω; brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1],G[i,ω] == 1);
@constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
@constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
#    @constraint(mp, tGnAnt[i in pData.II, ω in Ω], t[i] <= (disData[ω].H - ϵ) + G[i,ω]*M);
@constraint(mp, tGnAnt[i in pData.II, ω in Ω], t[i] <= disData[ω].H + G[i,ω]*M[i]);
@constraint(mp, tFnAnt[i in pData.II, ω in Ω], t[i] >= disData[ω].H - F[i,ω]*disData[ω].H);
@constraint(mp, FGnAnt[i in pData.II, ω in Ω], F[i,ω] + G[i,ω] == 1);

lb = -10000;
ub = 10000;

Mω = Dict();
for ω in Ω
    # Mω[ω] = 10*M;
    dDω = disData[ω];
    for i in pData.II
        Mω[i,ω] = lpSolveO(pData,dDω,i);
    end
end

TOL = 1e-6;

while abs((ub - lb)/lb) > 1e-6
    solve(mp);
    lbTemp = getobjectivevalue(mp);
    if lbTemp > lb
        lb = lbTemp;
    end
    xhat = getvalue(mp[:x]);
    that = getvalue(mp[:t]);
    Ghat = getvalue(mp[:G]);
    ubTemp = ubCal(pData,disData,Ω,xhat,that);
    if ubTemp < ub
        ub = ubTemp;
    end
    println(lb," ",ub);

    πr = Dict();
    λr = Dict();
    γr = Dict();
    vr = Dict();

    for ω in Ω
        sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0, CPX_PARAM_SIMDISPLAY = 0,CPX_PARAM_MIPDISPLAY = 0));
        @variable(sp, tω[i in pData.II] >= 0);
        @variable(sp, 0 <= xω[i in pData.II, j in pData.Ji[i]] <= 1);
        @variable(sp, 0 <= sω[i in pData.II, j in pData.Ji[i]] <= 1);

        @constraint(sp, durationConstr[k in pData.K], tω[k[2]] - tω[k[1]] + sum(pData.D[k[1]]*pData.eff[k[1]][j]*xω[k[1],j] for j in pData.Ji[k[1]])
            + sum(disData[ω].d[k[1]]*pData.eff[k[1]][j]*sω[k[1],j] for j in pData.Ji[k[1]])>= pData.D[k[1]] + disData[ω].d[k[1]]*getvalue(G[k[1],ω]));
        @constraint(sp, tGbound[i in pData.II], tω[i] >= disData[ω].H*Ghat[i,ω]);
        @constraint(sp, xConstr[i in pData.II], sum(xω[i,j] for j in pData.Ji[i]) <= 1);
        @constraint(sp, budgetConstr, sum(sum(xω[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
        @constraint(sp, tGnAnt1[i in pData.II], tω[i] <= that[i] + Ghat[i,ω]*Mω[i,ω]);
        @constraint(sp, tGnAnt2[i in pData.II], tω[i] >= that[i] - Ghat[i,ω]*M[i]);
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
                    disData[ω].H*getdual(tGbound)[i] + Mω[i,ω]*(getdual(tGnAnt1)[i]) - M[i]*(getdual(tGnAnt2)[i]) +
                    sum(getdual(xGnAnt1)[i,j] - getdual(xGnAnt2)[i,j] + getdual(sConstr2)[i,j] + getdual(sConstr3)[i,j] for j in pData.Ji[i]);
            else
                γr[ω][i] = disData[ω].H*getdual(tGbound)[i] + Mω[i,ω]*(getdual(tGnAnt1)[i]) - M[i]*(getdual(tGnAnt2)[i]);
            end
        end
        @constraint(mp, θ[ω] >= vω + sum(πr[ω][i]*(t[i] - that[i]) + γr[ω][i]*(G[i,ω] - Ghat[i,ω]) +
            sum(λr[ω][i,j]*(x[i,j] - xhat[i,j]) for j in pData.Ji[i]) for i in pData.II));
    end
    # @constraint(mp, θ >= sum(vr[ω]*disData[ω].prDis for ω in Ω) + sum(sum(πr[ω][i]*disData[ω].prDis for ω in Ω)*(t[i] - that[i]) + sum(γr[ω][i]*disData[ω].prDis*(G[i,ω] - Ghat[i,ω]) for ω in Ω) +
    #     sum(sum(λr[ω][i,j]*disData[ω].prDis for ω in Ω)*(x[i,j] - xhat[i,j]) for j in pData.Ji[i]) for i in pData.II));
end
