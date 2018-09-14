function genDisjunctive(pData,dDω,cutSetω,leafNodes,tm,xm,ts,xs,gs,ss,M,Mt)
    dp = Model(solver = GurobiSolver());
    # the number of disjunctive sections
    noT = length(leafNodes);
    nSet = 1:noT;
    gub = Dict();
    glb = Dict();
    for i in pData.II
        for n in nSet
            gub[i,n] = 1;
            glb[i,n] = 0;
            if i in leafNodes[n].zeroSet
                gub[i,n] = 0;
            end
            if i in leafNodes[n].oneSet
                glb[i,n] = 1;
            end
        end
    end
    # build the CGLP
    @variable(dp, -1 <= πt[i in pData.II] <= 1);
    @variable(dp, -1 <= πthat[i in pData.II] <= 1);
    @variable(dp, -1 <= λx[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(dp, -1 <= λxhat[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(dp, -1 <= γg[i in pData.II] <= 1);
    @variable(dp, -1 <= νs[i in pData.II, j in pData.Ji[i]] <= 1);

    @variable(dp, μsep[k in pData.K, n in nSet] >= 0);
    @variable(dp, μsephat[k in pData.K, n in nSet] >= 0);
    @variable(dp, μxlim[i in pData.II, n in nSet] <= 0);
    @variable(dp, μxhatlim[i in pData.II, n in nSet] <= 0);
    @variable(dp, μbudget[n in nSet] <= 0);
    @variable(dp, μbudgethat[n in nSet] <= 0);

    @variable(dp, μsg1[i in pData.II, j in pData.Ji[i], n in nSet] <= 0);
    @variable(dp, μsg2[i in pData.II, j in pData.Ji[i], n in nSet] <= 0);
    @variable(dp, μsg3[i in pData.II, j in pData.Ji[i], n in nSet] <= 0);
    @variable(dp, μgg[k in pData.K, n in nSet] <= 0);

    @variable(dp, μtGB[i in pData.II, n in nSet] <= 0);
    @variable(dp, μtG1[i in pData.II, n in nSet] >= 0);
    @variable(dp, μtG2[i in pData.II, n in nSet] <= 0);
    @variable(dp, μxG1[i in pData.II, j in pData.Ji[i], n in nSet] >= 0);
    @variable(dp, μxG2[i in pData.II, j in pData.Ji[i], n in nSet] <= 0);
    @variable(dp, μAnt1[i in pData.II, n in nSet] <= 0);
    @variable(dp, μAnt2[i in pData.II, n in nSet] <= 0);

    @variable(dp, μcut[l in 1:length(cutSetω), n in nSet] <= 0);

    @variable(dp, μtub[i in pData.II, n in nSet] <= 0);
    @variable(dp, μxub[i in pData.II, j in pData.Ji[i], n in nSet] <= 0);
    @variable(dp, μthatub[i in pData.II, n in nSet] <= 0);
    @variable(dp, μxhatub[i in pData.II, j in pData.Ji[i], n in nSet] <= 0);
    @variable(dp, μsub[i in pData.II, j in pData.Ji[i], n in nSet] <= 0);
    @variable(dp, μglb[i in pData.II, n in nSet] <= 0);
    @variable(dp, μgub[i in pData.II, n in nSet] <= 0);
    @variable(dp, μasum);

    @constraint(dp, tConstr[i in pData.II, n in nSet], -πt[i] + sum(μsep[k,n] for k in pData.K if k[1] == i) -
        sum(μsep[k,n] for k in pData.K if k[2] == i) + μtGB[i,n] - μtG1[i,n] - μtG2[i,n] - μtub[i,n] +
        sum(μcut[l,n]*cutSetω[l][3][i] for l in 1:length(cutSetω)) >= 0);
    @constraint(dp, xConstr[i in pData.II, j in pData.Ji[i], n in nSet], -λx[i,j] - sum(pData.D[i]*pData.eff[i][j]*μsep[k,n] for k in pData.K if k[1] == i) -
        μxlim[i,n] - μbudget[n] + μsg2[i,j,n] - μsg3[i,j,n] - μxG1[i,j,n] - μxG2[i,j,n] - μxub[i,j,n] +
        sum(μcut[l,n]*cutSetω[l][4][i,j] for l in 1:length(cutSetω)) >= 0);
    @constraint(dp, thatConstr[i in pData.II, n in nSet], -πthat[i] + sum(μsephat[k,n] for k in pData.K if k[1] == i) -
        sum(μsephat[k,n] for k in pData.K if k[2] == i) + μtG1[i,n] + μtG2[i,n] + μAnt1[i,n] - μAnt2[i,n] - μthatub[i,n] +
        sum(μcut[l,n]*cutSetω[l][1][i] for l in 1:length(cutSetω)) >= 0);
    @constraint(dp, xhatConstr[i in pData.II, j in pData.Ji[i], n in nSet], -λxhat[i,j] - sum(pData.D[i]*pData.eff[i][j]*μsephat[k,n] for k in pData.K if k[1] == i) -
        μxhatlim[i,n] - μbudgethat[n] + μxG1[i,j,n] + μxG2[i,j,n] - μxhatub[i,j,n] + sum(μcut[l,n]*cutSetω[l][2][i,j] for l in 1:length(cutSetω)) >= 0);
    @constraint(dp, gConstr[i in pData.II, n in nSet], -γg[i] + sum(dDω.d[i]*μsep[k,n] - μgg[k,n] for k in pData.K if k[1] == i) +
        sum(μsg1[i,j,n] - μsg3[i,j,n] for j in pData.Ji[i]) + sum(μgg[k,n] for k in pData.K if k[2] == i) - dDω.H*μtGB[i,n] -
        M*μtG1[i,n] + M*μtG2[i,n] + sum(-μxG1[i,j,n] + μxG2[i,j,n] for j in pData.Ji[i]) - Mt*μAnt1[i,n] + Mt*μAnt2[i,n] -
        μgub[i,n] + μglb[i,n] + sum(μcut[l,n]*cutSetω[l][5][i] for l in 1:length(cutSetω)) >= 0);
    @constraint(dp, sConstr[i in pData.II, j in pData.Ji[i], n in nSet], -νs[i,j] - sum(μsep[k,n]*dDω.d[i]*pData.eff[i][j] for k in pData.K if k[1] == i) -
        μsg1[i,j,n] - μsg2[i,j,n] + μsg3[i,j,n] - μsub[i,j,n] + sum(μcut[l,n]*cutSetω[l][6][i,j] for l in 1:length(cutSetω)) >= 0);
    @constraint(dp, aConstr[n in nSet], sum(pData.D[k[1]]*(μsep[k,n] + μsephat[k,n]) for k in pData.K) + pData.B*(μbudget[n] + μbudgethat[n]) +
        sum(μxlim[i,n] + μxhatlim[i,n] + (Mt - dDω.H)*μAnt1[i,n] + dDω.H*μAnt2[i,n] + M*μtub[i,n] + gub[i,n]*μgub[i,n] - glb[i,n]*μglb[i,n] +
        Mt*μthatub[i,n] for i in pData.II) - μasum - sum(μcut[l,n]*cutSetω[l][7] for l in 1:length(cutSetω)) +
        sum(sum(μsg3[i,j,n] + μsub[i,j,n] + μxub[i,j,n] + μxhatub[i,j,n] for j in pData.Ji[i]) for i in pData.II) >= 0);

    @objective(dp, Max, sum(πt[i]*ts[i] + πthat[i]*tm[i] + γg[i]*gs[i] +
        sum(λx[i,j]*xs[i,j] + λxhat[i,j]*xm[i,j] + νs[i,j]*ss[i,j] for j in pData.Ji[i]) for i in pData.II) + μasum);

    solve(dp);
    v = getobjectivevalue(dp);
    πDict = Dict();
    λDict = Dict();
    γDict = Dict();
    νDict = Dict();
    πhatDict = Dict();
    λhatDict = Dict();

    for i in pData.II
        πDict[i] = getvalue(dp[:πt][i]);
        πhatDict[i] = getvalue(dp[:πthat][i]);
        γDict[i] = getvalue(dp[:γg][i]);
        for j in pData.Ji[i]
            λDict[i,j] = getvalue(dp[:λx][i,j]);
            λhatDict[i,j] = getvalue(dp[:λxhat][i,j]);
            νDict[i,j] = getvalue(dp[:νs][i,j]);
        end
    end
    vo = getvalue(dp[:μasum]);
    # for i in pData.II
    #     vo -= πDict[i]*ts[i] + γDict[i]*gs[i] + πhatDict[i]*tm[i];
    #     for j in pData.Ji[i]
    #         vo -= λDict[i,j]*xs[i,j] + νDict[i,j]*ss[i,j] + λhatDict[i,j]*xm[i,j]
    #     end
    # end
    return vo,v,πDict,λDict,γDict,νDict,πhatDict,λhatDict;
end

function genDisjunctiveP(pData,dDω,cutSetω,leafNodes,tm,xm,ts,xs,gs,ss,M,Mt)
    noT = length(leafNodes);
    nSet = 1:noT;
    gub = Dict();
    glb = Dict();
    for i in pData.II
        for n in nSet
            gub[i,n] = 1;
            glb[i,n] = 0;
            if i in leafNodes[n].zeroSet
                gub[i,n] = 0;
            end
            if i in leafNodes[n].oneSet
                glb[i,n] = 1;
            end
        end
    end

    mp = Model(solver = GurobiSolver());

    @variable(mp, tplus[i in pData.II] >= 0);
    @variable(mp, tminus[i in pData.II] >= 0);
    @variable(mp, xplus[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(mp, xminus[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(mp, Gplus[i in pData.II] >= 0);
    @variable(mp, Gminus[i in pData.II] >= 0);
    @variable(mp, splus[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(mp, sminus[i in pData.II, j in pData.Ji[i]] >= 0);

    @variable(mp, thatplus[i in pData.II] >= 0);
    @variable(mp, thatminus[i in pData.II] >= 0);
    @variable(mp, xhatplus[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(mp, xhatminus[i in pData.II, j in pData.Ji[i]] >= 0);

    @variable(mp, a[n in nSet] >= 0);
    @variable(mp, t[i in pData.II, n in nSet] >= 0);
    @variable(mp, that[i in pData.II, n in nSet] >= 0);
    @variable(mp, x[i in pData.II, j in pData.Ji[i], n in nSet] >= 0);
    @variable(mp, xhat[i in pData.II, j in pData.Ji[i], n in nSet] >= 0);
    @variable(mp, G[i in pData.II, n in nSet] >= 0);
    @variable(mp, s[i in pData.II, j in pData.Ji[i], n in nSet] >= 0);

    @constraint(mp, tSum[i in pData.II], sum(t[i,n] for n in nSet) + tplus[i] - tminus[i] == ts[i]);
    @constraint(mp, xSum[i in pData.II, j in pData.Ji[i]], sum(x[i,j,n] for n in nSet) + xplus[i,j] - xminus[i,j] == xs[i,j]);
    @constraint(mp, thatSum[i in pData.II], sum(that[i,n] for n in nSet) + thatplus[i] - thatminus[i] == tm[i]);
    @constraint(mp, xhatSum[i in pData.II, j in pData.Ji[i]], sum(xhat[i,j,n] for n in nSet) + xhatplus[i,j] - xhatminus[i,j] == xm[i,j]);
    @constraint(mp, GSum[i in pData.II], sum(G[i,n] for n in nSet) + Gplus[i] - Gminus[i] == gs[i]);
    @constraint(mp, sSum[i in pData.II, j in pData.Ji[i]], sum(s[i,j,n] for n in nSet) + splus[i,j] - sminus[i,j] == ss[i,j]);
    @constraint(mp, aSum, sum(a[n] for n in nSet) == 1);

    @constraint(mp, sephat[k in pData.K, n in nSet], that[k[2],n] - that[k[1],n] >= pData.D[k[1]]*a[n] - sum(pData.eff[k[1]][j]*pData.D[k[1]]*xhat[k[1],j,n] for j in pData.Ji[k[1]]));
    @constraint(mp, xlimhat[i in pData.II, n in nSet], sum(xhat[i,j,n] for j in pData.Ji[i]) <= a[n]);
    @constraint(mp, budgethat[n in nSet], sum(sum(xhat[i,j,n] for j in pData.Ji[i]) for i in pData.II) <= pData.B*a[n]);

    @constraint(mp, sep[k in pData.K, n in nSet], t[k[2],n] - t[k[1],n] >= a[n]*pData.D[k[1]] + dDω.d[k[1]]*G[k[1],n] -
        sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j,n] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j,n] for j in pData.Ji[k[1]]));
    @constraint(mp, xlim[i in pData.II, n in nSet], sum(x[i,j,n] for j in pData.Ji[i]) <= a[n]);
    @constraint(mp, budget[n in nSet], sum(sum(x[i,j,n] for j in pData.Ji[i]) for i in pData.II) <= pData.B*a[n]);
    @constraint(mp, sGconstr1[i in pData.II, j in pData.Ji[i], n in nSet], s[i,j,n] <= G[i,n]);
    @constraint(mp, sGconstr2[i in pData.II, j in pData.Ji[i], n in nSet], s[i,j,n] <= x[i,j,n]);
    @constraint(mp, sGconstr3[i in pData.II, j in pData.Ji[i], n in nSet], s[i,j,n] >= G[i,n] + x[i,j,n] - a[n]);
    @constraint(mp, GGconstr[k in pData.K, n in nSet], G[k[1],n] <= G[k[2],n]);
    @constraint(mp, tGconstrB[i in pData.II, n in nSet], t[i,n] >= dDω.H*G[i,n]);
    @constraint(mp, tGconstr1[i in pData.II, n in nSet], t[i,n] + G[i,n]*M >= that[i,n]);
    @constraint(mp, tGconstr2[i in pData.II, n in nSet], t[i,n] - G[i,n]*M <= that[i,n]);
    @constraint(mp, xGconstr1[i in pData.II, j in pData.Ji[i], n in nSet], x[i,j,n] + G[i,n] >= xhat[i,j,n]);
    @constraint(mp, xGconstr2[i in pData.II, j in pData.Ji[i], n in nSet], x[i,j,n] - G[i,n] <= xhat[i,j,n]);
    @constraint(mp, tGAnt1[i in pData.II, n in nSet], Mt*G[i,n] - that[i,n] <= (Mt - dDω.H)*a[n]);
    @constraint(mp, tGAnt2[i in pData.II, n in nSet], -Mt*G[i,n] + that[i,n] <= dDω.H*a[n]);

    @constraint(mp, prevCuts[l in 1:length(cutSetω), n in nSet], a[n]*cutSetω[l][7] <=
        sum(cutSetω[l][1][i]*tm[i,n] + cutSetω[l][3][i]*t[i,n] + cutSetω[5][i]*G[i,n] +
        sum(cutSetω[l][2][i,j]*xm[i,j,n] + cutSetω[l][4][i,j]*x[i,j,n] + cutSetω[l][6]*s[i,j,n] for j in pData.Ji[i]) for i in pData.II));

    @constraint(mp, tub[i in pData.II, n in nSet], t[i,n] <= M*a[n]);
    @constraint(mp, thatub[i in pData.II, n in nSet], that[i,n] <= Mt*a[n]);
    @constraint(mp, xub[i in pData.II, j in pData.Ji[i], n in nSet], x[i,j,n] <= a[n]);
    @constraint(mp, xhatub[i in pData.II, j in pData.Ji[i], n in nSet], xhat[i,j,n] <= a[n]);
    @constraint(mp, sub[i in pData.II, j in pData.Ji[i], n in nSet], s[i,j,n] <= a[n]);
    @constraint(mp, Gub[i in pData.II, n in nSet], G[i,n] <= gub[i,n]*a[n]);
    @constraint(mp, Glb[i in pData.II, n in nSet], G[i,n] >= glb[i,n]*a[n]);

    @objective(mp, Min, sum(tplus[i] + tminus[i] + Gplus[i] + Gminus[i] + thatplus[i] + thatminus[i] +
        sum(xplus[i,j] + xminus[i,j] + xhatplus[i,j] + xhatminus[i,j] + splus[i,j] + sminus[i,j] for j in pData.Ji[i]) for i in pData.II));

    return mp;
end
