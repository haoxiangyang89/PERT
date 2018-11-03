# This is the extensive formulation of PERT (1st case)
# the disruption does not affect the activities that have not been started

function obtainExtSols(mext,ω)
    solve(mext);
    texts = Dict();
    xexts = Dict();
    gexts = Dict();
    sexts = Dict();
    for i in pData.II
        texts[i] = getvalue(mext[:t][i,ω]);
        gexts[i] = getvalue(mext[:G][i,ω]);
        for j in pData.Ji[i]
            xexts[i,j] = getvalue(mext[:x][i,j,ω]);
            sexts[i,j] = getvalue(mext[:s][i,j,ω]);
        end
    end
    return texts,xexts,gexts,sexts;
end

function extForm(pData,disData,Ω,TL = Inf)
    M = Dict();
    for ω in Ω
        #M[ω] = sum(max(pData.D[i],pData.D[i]+disData[ω].d[i]) for i in pData.II if i != 0);
        M[ω] = 300;
    end

    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-7, TimeLimit = TL));
    # mp = Model(solver = CplexSolver(CPX_PARAM_EPRHS = 1e-7,CPX_PARAM_EPINT = 1e-7,CPX_PARAM_TILIM = TL));
    @variable(mp,t0[i in pData.II] >= 0);
    @variable(mp,0 <= x0[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(mp,t[i in pData.II, ω in Ω] >= 0);
    @variable(mp,0 <= x[i in pData.II, j in pData.Ji[i], ω in Ω] <= 1);
    @variable(mp,G[i in pData.II, ω in Ω], Bin);
    @variable(mp,0 <= s[i in pData.II, j in pData.Ji[i], ω in Ω] <= 1);

    @constraint(mp, FConstr[i in pData.II, ω in Ω], disData[ω].H - (1-G[i,ω])*M[ω] <= t0[i]);
    @constraint(mp, GConstr[i in pData.II, ω in Ω], disData[ω].H - 1e-6 + G[i,ω]*M[ω] >= t0[i]);
    @constraint(mp, tConstr1[i in pData.II, ω in Ω], t[i,ω] + G[i,ω]*M[ω] >= t0[i]);
    @constraint(mp, tConstr2[i in pData.II, ω in Ω], t[i,ω] - G[i,ω]*M[ω] <= t0[i]);
    @constraint(mp, tConstr3[i in pData.II, ω in Ω], t[i,ω] >= disData[ω].H * G[i,ω]);
    @constraint(mp, xConstr1[i in pData.II, j in pData.Ji[i], ω in Ω], x[i,j,ω] + G[i,ω] >= x0[i,j]);
    @constraint(mp, xConstr2[i in pData.II, j in pData.Ji[i], ω in Ω], x[i,j,ω] - G[i,ω] <= x0[i,j]);
    @constraint(mp, durationConstr1[k in pData.K, ω in Ω], t[k[2],ω] - t[k[1],ω] >= pData.D[k[1]] + disData[ω].d[k[1]]*G[k[1],ω]
                      - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j,ω] + disData[ω].d[k[1]]*pData.eff[k[1]][j]*s[k[1],j,ω] for j in pData.Ji[k[1]]));
    @constraint(mp, durationConstr2[k in pData.K], t0[k[2]] - t0[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x0[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II,ω in Ω], sum(x[i,j,ω] for j in pData.Ji[i]) <= 1);
    @constraint(mp, GGcons1[i in pData.II, ω in 1:(length(Ω) - 1)], G[i,ω] >= G[i,ω + 1]);
    @constraint(mp, GGcons2[i in pData.K, ω in Ω], G[k[2],ω] >= G[k[1],ω]);
    @constraint(mp, budgetConstr[ω in Ω], sum(sum(pData.b[i][j]*x[i,j,ω] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr0[i in pData.II], sum(x0[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(mp, budgetConstr0, sum(sum(pData.b[i][j]*x0[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, Slinear1[i in pData.II, j in pData.Ji[i], ω in Ω], s[i,j,ω] <= G[i,ω]);
    @constraint(mp, Slinear2[i in pData.II, j in pData.Ji[i], ω in Ω], s[i,j,ω] <= x[i,j,ω]);
    @constraint(mp, Slinear3[i in pData.II, j in pData.Ji[i], ω in Ω], s[i,j,ω] >= x[i,j,ω] - 1 + G[i,ω]);

    @objective(mp,Min,pData.p0*t0[0] + sum(disData[ω].prDis*t[0,ω] for ω in Ω));

    solve(mp);
    text = getvalue(mp[:t0]);
    xext = getvalue(mp[:x0]);
    gext = getvalue(mp[:G]);
    fext = getobjectivevalue(mp);
    return text,xext,fext,gext,mp;
end

function extForm_cheat(pData,disData,Ω,prec = 1e-4,TL = Inf)
    M = Dict();
    for ω in Ω
        M[ω] = sum(max(pData.D[i],pData.D[i]+disData[ω].d[i]) for i in pData.II if i != 0);
    end

    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-7, TimeLimit = TL,MIPGap = prec));
    # mp = Model(solver = CplexSolver(CPX_PARAM_EPRHS = 1e-7,CPX_PARAM_EPINT = 1e-7,CPX_PARAM_TILIM = TL));
    @variable(mp,t0[i in pData.II] >= 0);
    @variable(mp,0 <= x0[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(mp,t[i in pData.II, ω in Ω] >= 0);
    @variable(mp,0 <= x[i in pData.II, j in pData.Ji[i], ω in Ω] <= 1);
    @variable(mp,G[i in pData.II, ω in Ω], Bin);
    @variable(mp,0 <= s[i in pData.II, j in pData.Ji[i], ω in Ω] <= 1);

    @constraint(mp, FConstr[i in pData.II, ω in Ω], disData[ω].H - (1 - G[i,ω])*M[ω] <= t0[i]);
    @constraint(mp, GConstr[i in pData.II, ω in Ω], disData[ω].H + G[i,ω]*M[ω] >= t0[i]);
    @constraint(mp, tConstr1[i in pData.II, ω in Ω], t[i,ω] + G[i,ω]*M[ω] >= t0[i]);
    @constraint(mp, tConstr2[i in pData.II, ω in Ω], t[i,ω] - G[i,ω]*M[ω] <= t0[i]);
    @constraint(mp, tConstr3[i in pData.II, ω in Ω], t[i,ω] >= disData[ω].H * G[i,ω]);
    @constraint(mp, xConstr1[i in pData.II, j in pData.Ji[i], ω in Ω], x[i,j,ω] + G[i,ω] >= x0[i,j]);
    @constraint(mp, xConstr2[i in pData.II, j in pData.Ji[i], ω in Ω], x[i,j,ω] - G[i,ω] <= x0[i,j]);
    @constraint(mp, durationConstr1[k in pData.K, ω in Ω], t[k[2],ω] - t[k[1],ω] >= pData.D[k[1]] + disData[ω].d[k[1]]*G[k[1],ω]
                      - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j,ω] + disData[ω].d[k[1]]*pData.eff[k[1]][j]*s[k[1],j,ω] for j in pData.Ji[k[1]]));
    @constraint(mp, durationConstr2[k in pData.K], t0[k[2]] - t0[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x0[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II,ω in Ω], sum(x[i,j,ω] for j in pData.Ji[i]) <= 1);
    @constraint(mp, GGcons1[i in pData.II, ω in 1:(length(Ω) - 1)], G[i,ω] >= G[i,ω + 1]);
    @constraint(mp, GGcons2[k in pData.K, ω in Ω], G[k[2],ω] >= G[k[1],ω]);
    @constraint(mp, budgetConstr[ω in Ω], sum(sum(pData.b[i][j]*x[i,j,ω] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr0[i in pData.II], sum(x0[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(mp, budgetConstr0, sum(sum(pData.b[i][j]*x0[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, Slinear1[i in pData.II, j in pData.Ji[i], ω in Ω], s[i,j,ω] <= G[i,ω]);
    @constraint(mp, Slinear2[i in pData.II, j in pData.Ji[i], ω in Ω], s[i,j,ω] <= x[i,j,ω]);
    @constraint(mp, Slinear3[i in pData.II, j in pData.Ji[i], ω in Ω], s[i,j,ω] >= x[i,j,ω] - 1 + G[i,ω]);

    # for the same H, bind the G variables
    HDiff = Dict();
    for ω in Ω
        if disData[ω].H in keys(HDiff)
            push!(HDiff[disData[ω].H], ω);
        else
            HDiff[disData[ω].H] = [ω];
        end
    end
    for hIter in keys(HDiff)
        if length(HDiff[hIter]) > 1
            for item in 2:length(HDiff[hIter])
                for i in pData.II
                    @constraint(mp, G[i,HDiff[hIter][item]] == G[i,HDiff[hIter][1]]);
                end
            end
        end
    end

    @objective(mp,Min,pData.p0*t0[0] + sum(disData[ω].prDis*t[0,ω] for ω in Ω));

    solve(mp);

    text = Dict();
    xext = Dict();
    gext = Dict();
    for i in pData.II
        text[i] = getvalue(mp[:t0][i]);
        for j in pData.Ji[i]
            xext[i,j] = getvalue(mp[:x0][i,j]);
        end
        for ω in Ω
            gext[i,ω] = getvalue(mp[:G][i,ω]);
        end
    end
    fext = getobjectivevalue(mp);
    return text,xext,fext,gext,mp;
end
