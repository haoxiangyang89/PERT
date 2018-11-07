# This is the script with all cut generation
function subRelax(pData,dDω,xhat,that,brInfo)
    M = sum(max(pData.D[i],pData.D[i]+dDω.d[i]) for i in pData.II if i != 0);

    # sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= F[i in pData.II] <= 1);
    @variable(sp, 0 <= G[i in pData.II] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0],dDω.H - (1 - G[i])*M <= that[i]);
    @constraint(sp, GCons[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0],dDω.H + G[i]*M >= that[i]);
    @constraint(sp, FFixed[i in pData.II; brInfo[findin(pData.II,i)[1]] == -1],G[i] == 0);
    @constraint(sp, GFixed[i in pData.II; brInfo[findin(pData.II,i)[1]] == 1],G[i] == 1);

    # add the predecessors and the successors logic constraints
    @constraint(sp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(sp, Min, t[0]);

    # obtain the dual variables for cuts
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v - λ(x - xhat) - π(t - that)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    for i in pData.II
        if brInfo[findin(pData.II,i)[1]] == 0
            πdict[i] = (getdual(sp[:FCons])[i] + getdual(sp[:GCons])[i] + getdual(sp[:tFnAnt1])[i] + getdual(sp[:tFnAnt2])[i]);
        else
            πdict[i] = (getdual(sp[:tFnAnt1])[i] + getdual(sp[:tFnAnt2])[i]);
        end
        for j in pData.Ji[i]
            λdict[i,j] = (getdual(sp[:xFnAnt1])[i,j] + getdual(sp[:xFnAnt2])[i,j]);
        end
    end
    cutGen = cutType(πdict,λdict,vk);
    return cutGen;
end

function subRelax_New(pData,dDω,xhat,that,brInfoω,M = Dict())
    # imbed the logic in the program
    if M == Dict()
        MM = iSolve_NC(pData,dDω,0,brInfoω);
        for i in pData.II
            M[i] = MM;
        end
    end

    # sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= G[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0],dDω.H - (1 - G[i])*M[i] <= that[i]);
    @constraint(sp, GCons[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0],dDω.H + G[i]*M[i] >= that[i]);

    # add the predecessors and the successors logic constraints
    @constraint(sp, GSuccessors[i in pData.II, k in pData.Succ[i]; (brInfoω[findin(pData.II,i)[1]] == 0)&(brInfoω[findin(pData.II,k)[1]] == 0)], G[i] <= G[k]);

    # add the basic sub problem constraints for the undecided activities
    @constraint(sp, tGbound1[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0],t[i] + G[i]*M[i] >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0],t[i] - G[i]*M[i] <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0],x[i,j] - G[i] <= xhat[i,j]);

    # add the logic constraints for the activities that have already started
    @constraint(sp, tFnAnt3[i in pData.II; brInfoω[findin(pData.II,i)[1]] == -1], t[i] == that[i]);
    @constraint(sp, xFnAnt3[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == -1],x[i,j] == xhat[i,j]);
    @constraint(sp, tGbound2[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 1],t[i] >= dDω.H);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K; brInfoω[findin(pData.II,k[1])[1]] == 0], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));
    @constraint(sp, durationConstr2[k in pData.K; brInfoω[findin(pData.II,k[1])[1]] == 1], t[k[2]] - t[k[1]] >= (pData.D[k[1]] + dDω.d[k[1]])*(1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(sp, durationConstr3[k in pData.K; brInfoω[findin(pData.II,k[1])[1]] == -1], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));

    @objective(sp, Min, t[0]);

    # obtain the dual variables for cuts
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v - λ(x - xhat) - π(t - that)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    for i in pData.II
        if brInfoω[findin(pData.II,i)[1]] == 0
            πdict[i] = (getdual(sp[:FCons])[i] + getdual(sp[:GCons])[i] + getdual(sp[:tFnAnt1])[i] + getdual(sp[:tFnAnt2])[i]);
            for j in pData.Ji[i]
                λdict[i,j] = (getdual(sp[:xFnAnt1])[i,j] + getdual(sp[:xFnAnt2])[i,j]);
            end
        elseif brInfoω[findin(pData.II,i)[1]] == -1
            πdict[i] = (getdual(sp[:tFnAnt3])[i]);
            for j in pData.Ji[i]
                λdict[i,j] = (getdual(sp[:xFnAnt3])[i,j]);
            end
        else
            πdict[i] = 0;
            for j in pData.Ji[i]
                λdict[i,j] = 0;
            end
        end
    end
    cutGen = cutType(πdict,λdict,vk);
    return cutGen;
end

function sub_Dual(pData,dDω,xhat,that,brInfoω,M = Dict())
    # imbed the logic in the program
    if M == Dict()
        MM = iSolve_NC(pData,dDω,0,brInfoω);
        for i in pData.II
            M[i] = MM;
        end
    end

    dp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));

    @variable(dp,λt0[k in pData.K; brInfoω[findin(pData.II,k[1])[1]] == 0] >= 0);
    @variable(dp,λth0[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0] >= 0);
    @variable(dp,λG1[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0] <= 0);
    @variable(dp,λG2[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0] >= 0);
    @variable(dp,λtM1[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0] <= 0);
    @variable(dp,λtM2[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0] >= 0);
    @variable(dp,λxM1[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0] <= 0);
    @variable(dp,λxM2[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0] >= 0);
    @variable(dp,λS1[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0] >= 0);
    @variable(dp,λS2[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0] <= 0);
    @variable(dp,λS3[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0] <= 0);
    @variable(dp,λGu[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0] <= 0);

    @variable(dp,λth1[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 1] >= 0);
    @variable(dp,λt1[k in pData.K; brInfoω[findin(pData.II,k[1])[1]] == 1] >= 0);

    @variable(dp,λtMe[i in pData.II; brInfoω[findin(pData.II,i)[1]] == -1]);
    @variable(dp,λxMe[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == -1]);
    @variable(dp,λte[k in pData.K; brInfoω[findin(pData.II,k[1])[1]] == -1] >= 0);

    @variable(dp,λB <= 0);
    @variable(dp,λxu[i in pData.II, j in pData.Ji[i]] <= 0);
    cDict = Dict();
    for i in pData.II
        if i == 0
            cDict[i] = 1;
        else
            cDict[i] = 0;
        end
    end

    @objective(dp, Max, sum(pData.D[k[1]]*λt0[k] for k in pData.K if brInfoω[findin(pData.II,k[1])[1]] == 0) +
        sum((that[i] + M[i] - dDω.H)*λG1[i] + (that[i] - dDω.H)*λG2[i] + that[i]*λtM1[i] + that[i]*λtM2[i] + λGu[i] for i in pData.II if brInfoω[findin(pData.II,i)[1]] == 0) +
        sum(sum(λxM1[i,j]*xhat[i,j] + λxM2[i,j]*xhat[i,j] - λS1[i,j] for j in pData.Ji[i]) for i in pData.II if brInfoω[findin(pData.II,i)[1]] == 0) +
        sum(dDω.H*λth1[i] for i in pData.II if brInfoω[findin(pData.II,i)[1]] == 1) + sum((pData.D[k[1]]+dDω.d[k[1]])*λt1[k] for k in pData.K if brInfoω[findin(pData.II,k[1])[1]] == 1) +
        sum(that[i]*λtMe[i] + sum(xhat[i,j]*λxMe[i,j] for j in pData.Ji[i]) for i in pData.II if brInfoω[findin(pData.II,i)[1]] == -1) +
        sum(pData.D[k[1]]*λte[k] for k in pData.K if brInfoω[findin(pData.II,k[1])[1]] == -1) +
        sum(sum(λxu[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λB);
    @constraint(dp,constrt0[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0],sum(-λt0[k] for k in pData.K if k[1] == i) + sum(λt0[k] for k in pData.K if (k[2] == i)&&(brInfoω[findin(pData.II,k[1])[1]] == 0)) +
        sum(λt1[k] for k in pData.K if (k[2] == i)&&(brInfoω[findin(pData.II,k[1])[1]] == 1)) + sum(λte[k] for k in pData.K if (k[2] == i)&&(brInfoω[findin(pData.II,k[1])[1]] == -1)) +
        λth0[i] + λtM1[i] + λtM2[i] <= cDict[i]);
    @constraint(dp,constrt1[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 1], sum(-λt1[k] for k in pData.K if k[1] == i) + sum(λt0[k] for k in pData.K if (k[2] == i)&&(brInfoω[findin(pData.II,k[1])[1]] == 0)) +
        sum(λt1[k] for k in pData.K if (k[2] == i)&&(brInfoω[findin(pData.II,k[1])[1]] == 1)) + sum(λte[k] for k in pData.K if (k[2] == i)&&(brInfoω[findin(pData.II,k[1])[1]] == -1)) +
        λth1[i] <= cDict[i]);
    @constraint(dp,constrte[i in pData.II; brInfoω[findin(pData.II,i)[1]] == -1],sum(-λte[k] for k in pData.K if k[1] == i) + sum(λt0[k] for k in pData.K if (k[2] == i)&&(brInfoω[findin(pData.II,k[1])[1]] == 0)) +
        sum(λt1[k] for k in pData.K if (k[2] == i)&&(brInfoω[findin(pData.II,k[1])[1]] == 1)) + sum(λte[k] for k in pData.K if (k[2] == i)&&(brInfoω[findin(pData.II,k[1])[1]] == -1)) +
        λtMe[i] <= cDict[i]);
    @constraint(dp,constrx0[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0], sum(pData.D[i]*pData.eff[i][j]*λt0[k] for k in pData.K if k[1] == i) +
        λxM1[i,j] + λxM2[i,j] - λS1[i,j] - λS3[i,j] + λB + λxu[i,j] <= 0);
    @constraint(dp,constrx1[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 1], sum((pData.D[i] + dDω.d[i])*pData.eff[i][j]*λt1[k] for k in pData.K if k[1] == i) +
        λB + λxu[i,j] <= 0);
    @constraint(dp,constrxe[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == -1],sum(pData.D[i]*pData.eff[i][j]*λte[k] for k in pData.K if k[1] == i) +
        λB + λxu[i,j] + λxMe[i,j] <= 0);
    @constraint(dp,constrG[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0], sum(-dDω.d[i]*λt0[k] for k in pData.K if k[1] == i) - dDω.H*λth0[i] + M[i]*λG1[i] + M[i]*λG2[i] -
        M[i]*λtM1[i] + M[i]*λtM2[i] - sum(λxM1[i,j] - λxM2[i,j] + λS1[i,j] + λS2[i,j] for j in pData.Ji[i]) + λGu[i] <= 0);
    @constraint(dp,constrS[i in pData.II, j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0], λS1[i,j] + λS2[i,j] + λS3[i,j] + sum(dDω.d[i]*pData.eff[i][j]*λt0[k] for k in pData.K if k[1] == i) <= 0);

    # obtain the dual variables for cuts
    solve(dp);
    vk = getobjectivevalue(dp);
    # the cut generated is θ >= v - λ(x - xhat) - π(t - that)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    for i in pData.II
        if brInfoω[findin(pData.II,i)[1]] == 0
            πdict[i] = (getvalue(dp[:λG1][i]) + getvalue(dp[:λG2][i]) + getvalue(dp[:λtM1][i]) + getvalue(dp[:λtM2][i]));
            for j in pData.Ji[i]
                λdict[i,j] = (getvalue(dp[:λxM1][i,j]) + getvalue(dp[:λxM2][i,j]));
            end
        elseif brInfoω[findin(pData.II,i)[1]] == -1
            πdict[i] = (getvalue(dp[:λtMe][i]));
            for j in pData.Ji[i]
                λdict[i,j] = (getvalue(dp[:λxMe][i,j]));
            end
        else
            πdict[i] = 0;
            for j in pData.Ji[i]
                λdict[i,j] = 0;
            end
        end
    end
    cutGen = cutType(πdict,λdict,vk);
    return cutGen;
end

# generate the
function subPull(pData,dDω,xhat,that,Ghatω,M = 9999999)

    # sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    #sp = Model(solver = CplexSolver(CPX_PARAM_EPINT = 1e-9,CPX_PARAM_EPRHS = 1e-9,CPX_PARAM_SCRIND = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0, IntFeasTol = 1e-9, FeasibilityTol = 1e-9));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    @variable(sp, s[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, 0 <= G[i in pData.II] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, Gcons[i in pData.II], G[i] == Ghatω[i]);
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(sp, Min, t[0]);

    # obtain the dual variables for cuts
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v + λ(x - xhat) + π(t - that) + γ(G - Ghat)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    γdict = Dict();             # dual for G
    for i in pData.II
        πdict[i] = (getdual(sp[:tFnAnt1])[i] + getdual(sp[:tFnAnt2])[i]);
        γdict[i] = getdual(sp[:Gcons][i]);
        for j in pData.Ji[i]
            λdict[i,j] = (getdual(sp[:xFnAnt1])[i,j] + getdual(sp[:xFnAnt2])[i,j]);
        end
    end
    return πdict,λdict,γdict,vk;
end

function subMixed(pData,dDω,xhat,that,Ghatω,ωCurr,Ω,M = 9999999,returnOpt = 0)
    sp = Model(solver = CplexSolver(CPX_PARAM_EPINT = 1e-9,CPX_PARAM_EPRHS = 1e-9,CPX_PARAM_SCRIND = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    @variable(sp, s[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, 0 <= G[i in pData.II] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, tGcons1[i in pData.II],dDω.H + G[i]*M >= that[i]);
    @constraint(sp, tGcons2[i in pData.II],dDω.H - (1 - G[i])*M <= that[i]);
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    # with the given Ghat
    @constraint(sp, GConstr[i in pData.II, j in pData.Succ[i]], G[i] <= G[j]);
    @constraint(sp, GGcons1[(i,ω) in keys(Ghatω);ω == ωCurr], G[i] == Ghatω[i,ω]);
    @constraint(sp, GGcons2[(i,ω) in keys(Ghatω);ω < ωCurr], G[i] <= Ghatω[i,ω]);
    @constraint(sp, GGcons3[(i,ω) in keys(Ghatω);ω > ωCurr], G[i] >= Ghatω[i,ω]);

    @objective(sp, Min, t[0]);

    # obtain the dual variables from solving the sub
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v + λ(x - xhat) + π(t - that) + γ(G - Ghat)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    γdict = Dict();             # dual for G
    for i in pData.II
        πdict[i] = (getdual(sp[:tFnAnt1][i]) + getdual(sp[:tFnAnt2][i]) + getdual(sp[:tGcons1][i]) + getdual(sp[:tGcons2][i]));
        for j in pData.Ji[i]
            λdict[i,j] = (getdual(sp[:xFnAnt1])[i,j] + getdual(sp[:xFnAnt2])[i,j]);
        end
    end
    for (i,ω) in keys(Ghat)
        if ω == ωCurr
            γdict[i,ω] = getdual(sp[:GGcons1][(i,ω)]);
        elseif ω < ωCurr
            γdict[i,ω] = getdual(sp[:GGcons2][(i,ω)]);
        else
            γdict[i,ω] = getdual(sp[:GGcons3][(i,ω)]);
        end
    end
    if returnOpt == 0
        return πdict,λdict,γdict,vk;
    else
        return πdict,λdict,γdict,vk,sp;
    end
end

function subF(pData,dDω,that,xhat,Fhat,M = 999999)
    sp = Model(solver = CplexSolver(CPX_PARAM_EPINT = 1e-9,CPX_PARAM_EPRHS = 1e-9,CPX_PARAM_SCRIND = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    @variable(sp, s[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, 0 <= G[i in pData.II] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, GFcons[i in pData.II, ω in Ω], G[i,ω] == sum(Fhat[i,ω1] for ω1 in 0:length(Ω) if ω1 >= ω));
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(sp, Min, t[0]);

    # obtain the dual variables for cuts
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v + λ(x - xhat) + π(t - that) + γ(G - Ghat)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    ηdict = Dict();             # dual for G
    for i in pData.II
        πdict[i] = (getdual(sp[:tFnAnt1])[i] + getdual(sp[:tFnAnt2])[i]);
        for ω in Ω
            ηdict[i] = getdual(sp[:GFcons][i,ω]);
        end
        for j in pData.Ji[i]
            λdict[i,j] = (getdual(sp[:xFnAnt1])[i,j] + getdual(sp[:xFnAnt2])[i,j]);
        end
    end
    return πdict,λdict,ηdict,vk;
end

function subLag(pData,dDω,xhat,that,Ghatω,γhat,M = 9999999)

    # sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    #sp = Model(solver = CplexSolver(CPX_PARAM_EPINT = 1e-9,CPX_PARAM_EPRHS = 1e-9,CPX_PARAM_SCRIND = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0, IntFeasTol = 1e-9, FeasibilityTol = 1e-9));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    @variable(sp, s[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, G[i in pData.II], Bin);

    # add the basic sub problem constraints
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(sp, Min, t[0] + sum(γhat[i]*(G[i] - Ghatω[i]) for i in pData.II));

    # obtain the dual variables for cuts
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v + λ(x - xhat) + π(t - that) + γ(G - Ghat)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    for i in pData.II
        πdict[i] = (getdual(sp[:tFnAnt1])[i] + getdual(sp[:tFnAnt2])[i]);
        for j in pData.Ji[i]
            λdict[i,j] = (getdual(sp[:xFnAnt1])[i,j] + getdual(sp[:xFnAnt2])[i,j]);
        end
    end
    return πdict,λdict,vk;
end

function subTight(pData,dDω,xhat,that,ubInfo,lbInfo,returnOpt = 0)
    M = sum(max(pData.D[i],pData.D[i]+dDω.d[i]) for i in pData.II if i != 0);

    # sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= G[i in pData.II] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, GCons11[i in pData.II; (dDω.H <= ubInfo[i])&(dDω.H > lbInfo[i])],G[i] <= (that[i] - dDω.H)/(dDω.H - lbInfo[i]) + 1);
    @constraint(sp, GCons12[i in pData.II; (dDω.H < ubInfo[i])&(dDω.H >= lbInfo[i])],G[i] >= (that[i] - dDω.H)/(ubInfo[i] - dDω.H));
    @constraint(sp, GCons21[i in pData.II; (dDω.H <= ubInfo[i])&(dDω.H >= lbInfo[i])],G[i] <= (that[i] - dDω.H)/1000 + 1);
    @constraint(sp, GCons22[i in pData.II; (dDω.H <= ubInfo[i])&(dDω.H >= lbInfo[i])],G[i] >= (that[i] - dDω.H)/1000);
    @constraint(sp, GFixed0[i in pData.II; dDω.H > ubInfo[i]],G[i] == 0);
    @constraint(sp, GFixed1[i in pData.II; dDω.H < lbInfo[i]],G[i] == 1);

    # add the predecessors and the successors logic constraints
    @constraint(sp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(sp, Min, t[0]);

    # obtain the dual variables for cuts
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v - λ(x - xhat) - π(t - that)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    for i in pData.II
        πdict[i] = (getdual(sp[:tFnAnt1])[i] + getdual(sp[:tFnAnt2])[i]);
        if (dDω.H <= ubInfo[i])&(dDω.H >= lbInfo[i])
            if (dDω.H <= ubInfo[i])&(dDω.H > lbInfo[i])
                πdict[i] += (getdual(sp[:GCons11])[i])/(dDω.H - lbInfo[i]);
            end
            if (dDω.H < ubInfo[i])&(dDω.H >= lbInfo[i])
                πdict[i] += (getdual(sp[:GCons12])[i])/(ubInfo[i] - dDω.H);
            end
            πdict[i] += (getdual(sp[:GCons21])[i])/1000 + (getdual(sp[:GCons22])[i])/1000;
        end
        for j in pData.Ji[i]
            λdict[i,j] = (getdual(sp[:xFnAnt1])[i,j] + getdual(sp[:xFnAnt2])[i,j]);
        end
    end
    if returnOpt == 0
        return πdict,λdict,vk;
    else
        return πdict,λdict,vk,sp;
    end
end

function sub_div(pData,dDω,ωCurr,that,xhat,yhat,divSet,ubInfo,lbInfo,M1 = 999999,returnOpt = 0)
    M = sum(max(pData.D[i],pData.D[i]+dDω.d[i]) for i in pData.II if i != 0);

    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= G[i in pData.II] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, GCons11[i in pData.II],G[i] <= (that[i] - dDω.H)/M1 + 1);
    @constraint(sp, GCons12[i in pData.II],G[i] >= (that[i] - dDω.H)/M1);
    @constraint(sp, GCons21[i in pData.II],G[i] <= (that[i] - dDω.H)/(dDω.H - lbInfo[i]) + 1);
    @constraint(sp, GCons22[i in pData.II],G[i] >= (that[i] - dDω.H)/(ubInfo[i] - dDω.H));
    @constraint(sp, GFixed0[i in pData.II],G[i] >= sum(yhat[i,par] for par in 1:length(divSet[i]) if ωCurr < divSet[i][par].startH));
    @constraint(sp, GFixed1[i in pData.II],G[i] <= 1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH]));

    # add the predecessors and the successors logic constraints
    @constraint(sp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(sp, Min, t[0]);

    # obtain the dual variables for cuts
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v - λ(x - xhat) - π(t - that)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    γdict = Dict();             # dual for y
    for i in pData.II
        πdict[i] = (getdual(sp[:GCons11][i]) + getdual(sp[:GCons12][i]))/M1 +
            getdual(sp[:GCons21][i])/(dDω.H - lbInfo[i]) + getdual(sp[:GCons22][i])/(ubInfo[i] - dDω.H) +
            (getdual(sp[:tFnAnt1][i]) + getdual(sp[:tFnAnt2][i]));
        for j in pData.Ji[i]
            λdict[i,j] = (getdual(sp[:xFnAnt1][i,j]) + getdual(sp[:xFnAnt2][i,j]));
        end
        for par in 1:length(divSet[i])
            if ωCurr < divSet[i][par].startH
                γdict[i,par] = getdual(sp[:GFixed0][i]);
            elseif dDω.H >= H[divSet[i][par].endH]
                γdict[i,par] = -getdual(sp[:GFixed1][i]);
            else
                γdict[i,par] = 0;
            end
        end
    end
    Ghat = Dict();
    for i in pData.II
        Ghat[i] = getvalue(sp[:G][i]);
    end
    if returnOpt == 0
        return πdict,λdict,γdict,vk,Ghat;
    else
        return πdict,λdict,γdict,vk,sp;
    end
end

function sub_divT(pData,dDω,ωCurr,that,xhat,yhat,divSet,H,M,returnOpt = 0)
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= G[i in pData.II] <= 1);
    @variable(sp, 0 <= Gy[i in pData.II,par in 1:length(divSet[i])] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, GyRelax1[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] <= G[i]);
    @constraint(sp, GyRelax2[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] <= yhat[i,par]);
    @constraint(sp, GyRelax3[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] >= G[i] + yhat[i,par] - 1);

    @constraint(sp, GCons1[i in pData.II],dDω.H*G[i] - sum(H[divSet[i][par].endH]*Gy[i,par] for par in 1:length(divSet[i])) <= dDω.H - that[i]);
    @constraint(sp, GCons2[i in pData.II],sum(H[divSet[i][par].startH]*Gy[i,par] for par in 1:length(divSet[i])) - dDω.H*G[i] >= -that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i])));
    @constraint(sp, GFixed0[i in pData.II],G[i] >= sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH]));
    @constraint(sp, GFixed1[i in pData.II],G[i] <= 1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH]));

    # add the predecessors and the successors logic constraints
    @constraint(sp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M[i] >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M[i] <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(sp, Min, t[0]);

    # obtain the dual variables for cuts
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v - λ(x - xhat) - π(t - that)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    γdict = Dict();             # dual for y
    for i in pData.II
        πdict[i] = -getdual(sp[:GCons1][i]) - getdual(sp[:GCons2][i]) +
            (getdual(sp[:tFnAnt1][i]) + getdual(sp[:tFnAnt2][i]));
        for j in pData.Ji[i]
            λdict[i,j] = (getdual(sp[:xFnAnt1][i,j]) + getdual(sp[:xFnAnt2][i,j]));
        end
        for par in 1:length(divSet[i])
            γdict[i,par] = H[divSet[i][par].startH]*getdual(sp[:GCons2][i]) +
                getdual(sp[:GyRelax2][i,par]) + getdual(sp[:GyRelax3][i,par]);
            if dDω.H <= H[divSet[i][par].startH]
                γdict[i,par] += getdual(sp[:GFixed0][i]);
            elseif dDω.H >= H[divSet[i][par].endH]
                γdict[i,par] -= getdual(sp[:GFixed1][i]);
            end
        end
    end
    Ghat = Dict();
    for i in pData.II
        Ghat[i] = getvalue(sp[:G][i]);
    end
    if returnOpt == 0
        return πdict,λdict,γdict,vk,Ghat;
    else
        return πdict,λdict,γdict,vk,sp;
    end
end

function sub_divTDual(pData,dDω,ωCurr,that,xhat,yhat,divSet,H,M,returnOpt = 0)
    # solve the subproblem by dual formulation
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, λFG1[i in pData.II, par in 1:length(divSet[i])] <= 0);
    @variable(sp, λFG2[i in pData.II, par in 1:length(divSet[i])] <= 0);
    @variable(sp, λFG3[i in pData.II, par in 1:length(divSet[i])] >= 0);
    @variable(sp, λHG1[i in pData.II] <= 0);
    @variable(sp, λHG2[i in pData.II] >= 0);
    @variable(sp, λGy1[i in pData.II] >= 0);
    @variable(sp, λGy2[i in pData.II] <= 0);
    @variable(sp, λGG[k in pData.K] <= 0);
    @variable(sp, λtG1[i in pData.II] >= 0);
    @variable(sp, λtG2[i in pData.II] >= 0);
    @variable(sp, λtG3[i in pData.II] <= 0);
    @variable(sp, λxG1[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, λxG2[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG1[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG2[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG3[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, λbudget <= 0);
    @variable(sp, λxub[i in pData.II] <= 0);
    @variable(sp, λdur[k in pData.K] >= 0);

    tBool = Dict();
    for i in pData.II
        if i != 0
            tBool[i] = 0;
        else
            tBool[i] = 1;
        end
    end

    # t constraint
    @constraint(sp, tConstr[i in pData.II], tBool[i] - λtG1[i] - λtG2[i] - λtG3[i] +
        sum(λdur[k] for k in pData.K if k[1] == i) - sum(λdur[k] for k in pData.K if k[2] == i) >= 0);
    # x constraint
    @constraint(sp, xConstr[i in pData.II, j in pData.Ji[i]], -λxG1[i,j] - λxG2[i,j] +
        λsG2[i,j] + λsG3[i,j] - pData.b[i][j]*λbudget - λxub[i] - sum(pData.D[i]*pData.eff[i][j]*λdur[k] for k in pData.K if k[1] == i) >= 0);
    # G constraint
    @constraint(sp, gConstr[i in pData.II], sum(λFG1[i,par] + λFG3[i,par] for par in 1:length(divSet[i])) - dDω.H*(λHG1[i] - λHG2[i]) -
        λGy1[i] - λGy2[i] + sum(λGG[k] for k in pData.K if k[2] == i) - sum(λGG[k] for k in pData.K if k[1] == i) + λtG1[i]*dDω.H -
        λtG2[i]*M[i] + λtG3[i]*M[i] + sum(-λxG1[i,j] + λxG2[i,j] + λsG1[i,j] + λsG3[i,j] for j in pData.Ji[i]) +
        sum(dDω.d[i]*λdur[k] for k in pData.K if k[1] == i) >= 0);
    # Gy constraint
    @constraint(sp, FConstr[i in pData.II, par in 1:length(divSet[i])], -λFG1[i,par] - λFG2[i,par] - λFG3[i,par] +
        H[divSet[i][par].endH]*λHG1[i] - H[divSet[i][par].startH]*λHG2[i] >= 0);
    # s constraint
    @constraint(sp, sConstr[i in pData.II, j in pData.Ji[i]], -λsG1[i,j] - λsG2[i,j] - λsG3[i,j] -
        sum(dDω.d[i]*pData.eff[i][j]*λdur[k] for k in pData.K if k[1] == i) >= 0);

    # objective function
    @objective(sp, Max, sum(sum(yhat[i,par]*λFG2[i,par] + (yhat[i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - that[i]) + λHG2[i]*(-that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(that[i]*(λtG2[i] + λtG3[i]) + sum(xhat[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K));

    solve(sp);
    vk = getobjectivevalue(sp);

    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    γdict = Dict();             # dual for y
    for i in pData.II
        πdict[i] = -getvalue(sp[:λHG1][i]) - getvalue(sp[:λHG2][i]) + getvalue(sp[:λtG2][i]) + getvalue(sp[:λtG3][i]);
        for j in pData.Ji[i]
            λdict[i,j] = getvalue(sp[:λxG1][i,j]) + getvalue(sp[:λxG2][i,j]);
        end
        for par in 1:length(divSet[i])
            γdict[i,par] = H[divSet[i][par].startH]*getvalue(sp[:λHG2][i]) +
                getvalue(sp[:λFG2][i,par]) + getvalue(sp[:λFG3][i,par]);
            if dDω.H <= H[divSet[i][par].startH]
                γdict[i,par] += getvalue(sp[:λGy1][i]);
            elseif dDω.H >= H[divSet[i][par].endH]
                γdict[i,par] -= getvalue(sp[:λGy2][i]);
            end
        end
    end
    Ghat = Dict();
    for i in pData.II
        Ghat[i] = -getdual(sp[:gConstr][i]);
    end

    if returnOpt == 0
        return πdict,λdict,γdict,vk,Ghat;
    else
        return πdict,λdict,γdict,vk,sp;
    end
end

function sub_divTDualT(pData,dDω,ωCurr,that,xhat,yhat,divSet,H,M,tcore,xcore,ycore,returnOpt = 0)
    # smp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    smp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(smp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(smp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(smp, 0 <= G[i in pData.II] <= 1);
    @variable(smp, 0 <= Gy[i in pData.II,par in 1:length(divSet[i])] <= 1);
    @variable(smp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    # add the basic sub problem constraints
    @constraint(smp, GyRelax1[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] <= G[i]);
    @constraint(smp, GyRelax2[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] <= yhat[i,par]);
    @constraint(smp, GyRelax3[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] >= G[i] + yhat[i,par] - 1);

    @constraint(smp, GCons1[i in pData.II],dDω.H*G[i] - sum(H[divSet[i][par].endH]*Gy[i,par] for par in 1:length(divSet[i])) <= dDω.H - that[i]);
    @constraint(smp, GCons2[i in pData.II],sum(H[divSet[i][par].startH]*Gy[i,par] for par in 1:length(divSet[i])) - dDω.H*G[i] >= -that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i])));
    @constraint(smp, GFixed0[i in pData.II],G[i] >= sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH]));
    @constraint(smp, GFixed1[i in pData.II],G[i] <= 1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH]));

    # add the predecessors and the successors logic constraints
    @constraint(smp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(smp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(smp, tFnAnt1[i in pData.II],t[i] + G[i]*M[i] >= that[i]);
    @constraint(smp, tFnAnt2[i in pData.II],t[i] - G[i]*M[i] <= that[i]);
    @constraint(smp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(smp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(smp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(smp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(smp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(smp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(smp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(smp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(smp, Min, t[0]);
    smpStatus = solve(smp);
    vhat1 = getobjectivevalue(smp);
    Ghat = Dict();
    for i in pData.II
        Ghat[i] = getvalue(smp[:G][i]);
    end
    λdict1 = Dict();             # dual for x
    πdict1 = Dict();             # dual for t
    γdict1 = Dict();             # dual for y
    for i in pData.II
        πdict1[i] = -getdual(smp[:GCons1][i]) - getdual(smp[:GCons2][i]) +
            (getdual(smp[:tFnAnt1][i]) + getdual(smp[:tFnAnt2][i]));
        for j in pData.Ji[i]
            λdict1[i,j] = (getdual(smp[:xFnAnt1][i,j]) + getdual(smp[:xFnAnt2][i,j]));
        end
        for par in 1:length(divSet[i])
            γdict1[i,par] = H[divSet[i][par].startH]*getdual(smp[:GCons2][i]) +
                getdual(smp[:GyRelax2][i,par]) + getdual(smp[:GyRelax3][i,par]);
            if dDω.H <= H[divSet[i][par].startH]
                γdict1[i,par] += getdual(smp[:GFixed0][i]);
            elseif dDω.H >= H[divSet[i][par].endH]
                γdict1[i,par] -= getdual(smp[:GFixed1][i]);
            end
        end
    end

    # solve the subproblem by dual formulation
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, λFG1[i in pData.II, par in 1:length(divSet[i])] <= 0);
    @variable(sp, λFG2[i in pData.II, par in 1:length(divSet[i])] <= 0);
    @variable(sp, λFG3[i in pData.II, par in 1:length(divSet[i])] >= 0);
    @variable(sp, λHG1[i in pData.II] <= 0);
    @variable(sp, λHG2[i in pData.II] >= 0);
    @variable(sp, λGy1[i in pData.II] >= 0);
    @variable(sp, λGy2[i in pData.II] <= 0);
    @variable(sp, λGG[k in pData.K] <= 0);
    @variable(sp, λtG1[i in pData.II] >= 0);
    @variable(sp, λtG2[i in pData.II] >= 0);
    @variable(sp, λtG3[i in pData.II] <= 0);
    @variable(sp, λxG1[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, λxG2[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG1[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG2[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG3[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, λbudget <= 0);
    @variable(sp, λxub[i in pData.II] <= 0);
    @variable(sp, λdur[k in pData.K] >= 0);

    tBool = Dict();
    for i in pData.II
        if i != 0
            tBool[i] = 0;
        else
            tBool[i] = 1;
        end
    end

    # t constraint
    @constraint(sp, tConstr[i in pData.II], tBool[i] - λtG1[i] - λtG2[i] - λtG3[i] +
        sum(λdur[k] for k in pData.K if k[1] == i) - sum(λdur[k] for k in pData.K if k[2] == i) >= 0);
    # x constraint
    @constraint(sp, xConstr[i in pData.II, j in pData.Ji[i]], -λxG1[i,j] - λxG2[i,j] +
        λsG2[i,j] + λsG3[i,j] - pData.b[i][j]*λbudget - λxub[i] - sum(pData.D[i]*pData.eff[i][j]*λdur[k] for k in pData.K if k[1] == i) >= 0);
    # G constraint
    @constraint(sp, gConstr[i in pData.II], sum(λFG1[i,par] + λFG3[i,par] for par in 1:length(divSet[i])) - dDω.H*(λHG1[i] - λHG2[i]) -
        λGy1[i] - λGy2[i] + sum(λGG[k] for k in pData.K if k[2] == i) - sum(λGG[k] for k in pData.K if k[1] == i) + λtG1[i]*dDω.H -
        λtG2[i]*M[i] + λtG3[i]*M[i] + sum(-λxG1[i,j] + λxG2[i,j] + λsG1[i,j] + λsG3[i,j] for j in pData.Ji[i]) +
        sum(dDω.d[i]*λdur[k] for k in pData.K if k[1] == i) >= 0);
    # Gy constraint
    @constraint(sp, FConstr[i in pData.II, par in 1:length(divSet[i])], -λFG1[i,par] - λFG2[i,par] - λFG3[i,par] +
        H[divSet[i][par].endH]*λHG1[i] - H[divSet[i][par].startH]*λHG2[i] >= 0);
    # s constraint
    @constraint(sp, sConstr[i in pData.II, j in pData.Ji[i]], -λsG1[i,j] - λsG2[i,j] - λsG3[i,j] -
        sum(dDω.d[i]*pData.eff[i][j]*λdur[k] for k in pData.K if k[1] == i) >= 0);

    @expression(sp, corePoint, sum(sum(yhat[i,par]*λFG2[i,par] + (yhat[i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - that[i]) + λHG2[i]*(-that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(that[i]*(λtG2[i] + λtG3[i]) + sum(xhat[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K));

    # objective function of the binary feasible solution should be the same
    @constraint(sp, binaryTight, sum(sum(yhat[i,par]*λFG2[i,par] + (yhat[i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - that[i]) + λHG2[i]*(-that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(that[i]*(λtG2[i] + λtG3[i]) + sum(xhat[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K) >= (1- 1e-5)*vhat1);

    # optimize the fractional solution's objective
    @objective(sp, Max, sum(sum(ycore[i,par]*λFG2[i,par] + (ycore[i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - tcore[i]) + λHG2[i]*(-tcore[i] + sum(ycore[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(ycore[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(ycore[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(tcore[i]*(λtG2[i] + λtG3[i]) + sum(xcore[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K));

    spStatus = solve(sp);
    vhat2 = getvalue(corePoint);

    if spStatus == :Optimal
        λdict2 = Dict();             # dual for x
        πdict2 = Dict();             # dual for t
        γdict2 = Dict();             # dual for y
        for i in pData.II
            πdict2[i] = -getvalue(sp[:λHG1][i]) - getvalue(sp[:λHG2][i]) + getvalue(sp[:λtG2][i]) + getvalue(sp[:λtG3][i]);
            for j in pData.Ji[i]
                λdict2[i,j] = getvalue(sp[:λxG1][i,j]) + getvalue(sp[:λxG2][i,j]);
            end
            for par in 1:length(divSet[i])
                γdict2[i,par] = H[divSet[i][par].startH]*getvalue(sp[:λHG2][i]) +
                    getvalue(sp[:λFG2][i,par]) + getvalue(sp[:λFG3][i,par]);
                if dDω.H <= H[divSet[i][par].startH]
                    γdict2[i,par] += getvalue(sp[:λGy1][i]);
                elseif dDω.H >= H[divSet[i][par].endH]
                    γdict2[i,par] -= getvalue(sp[:λGy2][i]);
                end
            end
        end
    else
        print("Error on $(ωCurr)"," ",spStatus);
    end

    if returnOpt == 0
        return πdict1,λdict1,γdict1,vhat1,Ghat,πdict2,λdict2,γdict2,vhat2;
    else
        return πdict1,λdict1,γdict1,vhat1,sp,πdict2,λdict2,γdict2,vhat2;
    end
end

function sub_divTDualT2(pData,dDω,ωCurr,that,xhat,yhat,divSet,H,M,tcore,xcore,ycore,returnOpt = 0)
    # Magnanti-Wong with a small perturbation
    # smp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    smp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(smp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(smp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(smp, 0 <= G[i in pData.II] <= 1);
    @variable(smp, 0 <= Gy[i in pData.II,par in 1:length(divSet[i])] <= 1);
    @variable(smp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    # add the basic sub problem constraints
    @constraint(smp, GyRelax1[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] <= G[i]);
    @constraint(smp, GyRelax2[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] <= yhat[i,par]);
    @constraint(smp, GyRelax3[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] >= G[i] + yhat[i,par] - 1);

    @constraint(smp, GCons1[i in pData.II],dDω.H*G[i] - sum(H[divSet[i][par].endH]*Gy[i,par] for par in 1:length(divSet[i])) <= dDω.H - that[i]);
    @constraint(smp, GCons2[i in pData.II],sum(H[divSet[i][par].startH]*Gy[i,par] for par in 1:length(divSet[i])) - dDω.H*G[i] >= -that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i])));
    @constraint(smp, GFixed0[i in pData.II],G[i] >= sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH]));
    @constraint(smp, GFixed1[i in pData.II],G[i] <= 1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH]));

    # add the predecessors and the successors logic constraints
    @constraint(smp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(smp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(smp, tFnAnt1[i in pData.II],t[i] + G[i]*M[i] >= that[i]);
    @constraint(smp, tFnAnt2[i in pData.II],t[i] - G[i]*M[i] <= that[i]);
    @constraint(smp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(smp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(smp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(smp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(smp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(smp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(smp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(smp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(smp, Min, t[0]);
    smpStatus = solve(smp);
    vhat = getobjectivevalue(smp);
    Ghat = Dict();
    for i in pData.II
        Ghat[i] = getvalue(smp[:G][i]);
    end

    # solve the subproblem by dual formulation
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, λFG1[i in pData.II, par in 1:length(divSet[i])] <= 0);
    @variable(sp, λFG2[i in pData.II, par in 1:length(divSet[i])] <= 0);
    @variable(sp, λFG3[i in pData.II, par in 1:length(divSet[i])] >= 0);
    @variable(sp, λHG1[i in pData.II] <= 0);
    @variable(sp, λHG2[i in pData.II] >= 0);
    @variable(sp, λGy1[i in pData.II] >= 0);
    @variable(sp, λGy2[i in pData.II] <= 0);
    @variable(sp, λGG[k in pData.K] <= 0);
    @variable(sp, λtG1[i in pData.II] >= 0);
    @variable(sp, λtG2[i in pData.II] >= 0);
    @variable(sp, λtG3[i in pData.II] <= 0);
    @variable(sp, λxG1[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, λxG2[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG1[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG2[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG3[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, λbudget <= 0);
    @variable(sp, λxub[i in pData.II] <= 0);
    @variable(sp, λdur[k in pData.K] >= 0);

    tBool = Dict();
    for i in pData.II
        if i != 0
            tBool[i] = 0;
        else
            tBool[i] = 1;
        end
    end

    # t constraint
    @constraint(sp, tConstr[i in pData.II], tBool[i] - λtG1[i] - λtG2[i] - λtG3[i] +
        sum(λdur[k] for k in pData.K if k[1] == i) - sum(λdur[k] for k in pData.K if k[2] == i) >= 0);
    # x constraint
    @constraint(sp, xConstr[i in pData.II, j in pData.Ji[i]], -λxG1[i,j] - λxG2[i,j] +
        λsG2[i,j] + λsG3[i,j] - pData.b[i][j]*λbudget - λxub[i] - sum(pData.D[i]*pData.eff[i][j]*λdur[k] for k in pData.K if k[1] == i) >= 0);
    # G constraint
    @constraint(sp, gConstr[i in pData.II], sum(λFG1[i,par] + λFG3[i,par] for par in 1:length(divSet[i])) - dDω.H*(λHG1[i] - λHG2[i]) -
        λGy1[i] - λGy2[i] + sum(λGG[k] for k in pData.K if k[2] == i) - sum(λGG[k] for k in pData.K if k[1] == i) + λtG1[i]*dDω.H -
        λtG2[i]*M[i] + λtG3[i]*M[i] + sum(-λxG1[i,j] + λxG2[i,j] + λsG1[i,j] + λsG3[i,j] for j in pData.Ji[i]) +
        sum(dDω.d[i]*λdur[k] for k in pData.K if k[1] == i) >= 0);
    # Gy constraint
    @constraint(sp, FConstr[i in pData.II, par in 1:length(divSet[i])], -λFG1[i,par] - λFG2[i,par] - λFG3[i,par] +
        H[divSet[i][par].endH]*λHG1[i] - H[divSet[i][par].startH]*λHG2[i] >= 0);
    # s constraint
    @constraint(sp, sConstr[i in pData.II, j in pData.Ji[i]], -λsG1[i,j] - λsG2[i,j] - λsG3[i,j] -
        sum(dDω.d[i]*pData.eff[i][j]*λdur[k] for k in pData.K if k[1] == i) >= 0);

    @expression(sp, corePoint, sum(sum(yhat[i,par]*λFG2[i,par] + (yhat[i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - that[i]) + λHG2[i]*(-that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(that[i]*(λtG2[i] + λtG3[i]) + sum(xhat[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K));

    # objective function of the binary feasible solution should be the same
    @constraint(sp, binaryTight, sum(sum(yhat[i,par]*λFG2[i,par] + (yhat[i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - that[i]) + λHG2[i]*(-that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(that[i]*(λtG2[i] + λtG3[i]) + sum(xhat[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K) >= (1- 1e-5)*vhat);

    # optimize the fractional solution's objective
    @objective(sp, Max, sum(sum(ycore[i,par]*λFG2[i,par] + (ycore[i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - tcore[i]) + λHG2[i]*(-tcore[i] + sum(ycore[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(ycore[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(ycore[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(tcore[i]*(λtG2[i] + λtG3[i]) + sum(xcore[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K));

    spStatus = solve(sp);

    if spStatus == :Optimal
        λdict = Dict();             # dual for x
        πdict = Dict();             # dual for t
        γdict = Dict();             # dual for y
        vhat = getvalue(corePoint);
        for i in pData.II
            πdict[i] = -getvalue(sp[:λHG1][i]) - getvalue(sp[:λHG2][i]) + getvalue(sp[:λtG2][i]) + getvalue(sp[:λtG3][i]);
            for j in pData.Ji[i]
                λdict[i,j] = getvalue(sp[:λxG1][i,j]) + getvalue(sp[:λxG2][i,j]);
            end
            for par in 1:length(divSet[i])
                γdict[i,par] = H[divSet[i][par].startH]*getvalue(sp[:λHG2][i]) +
                    getvalue(sp[:λFG2][i,par]) + getvalue(sp[:λFG3][i,par]);
                if dDω.H <= H[divSet[i][par].startH]
                    γdict[i,par] += getvalue(sp[:λGy1][i]);
                elseif dDω.H >= H[divSet[i][par].endH]
                    γdict[i,par] -= getvalue(sp[:λGy2][i]);
                end
            end
        end
    else
        vhat = getobjectivevalue(smp);
        λdict = Dict();             # dual for x
        πdict = Dict();             # dual for t
        γdict = Dict();             # dual for y
        for i in pData.II
            πdict[i] = -getdual(smp[:GCons1][i]) - getdual(smp[:GCons2][i]) +
                (getdual(smp[:tFnAnt1][i]) + getdual(smp[:tFnAnt2][i]));
            for j in pData.Ji[i]
                λdict[i,j] = (getdual(smp[:xFnAnt1][i,j]) + getdual(smp[:xFnAnt2][i,j]));
            end
            for par in 1:length(divSet[i])
                γdict[i,par] = H[divSet[i][par].startH]*getdual(smp[:GCons2][i]) +
                    getdual(smp[:GyRelax2][i,par]) + getdual(smp[:GyRelax3][i,par]);
                if dDω.H <= H[divSet[i][par].startH]
                    γdict[i,par] += getdual(smp[:GFixed0][i]);
                elseif dDω.H >= H[divSet[i][par].endH]
                    γdict[i,par] -= getdual(smp[:GFixed1][i]);
                end
            end
        end
    end

    if returnOpt == 0
        return πdict,λdict,γdict,vhat,Ghat;
    else
        return πdict,λdict,γdict,vhat,sp;
    end
end

function sub_divTDualT3(pData,dDω,ωCurr,that,xhat,yhat,divSet,H,M,tcoreList,xcoreList,ycoreList,returnOpt = 0)
    # give a list of previous recorded mip feasible solution
    # optimize over a weighted sum of the cut value at those previous solutions
    smp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    @variable(smp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(smp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(smp, 0 <= G[i in pData.II] <= 1);
    @variable(smp, 0 <= Gy[i in pData.II,par in 1:length(divSet[i])] <= 1);
    @variable(smp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    # add the basic sub problem constraints
    @constraint(smp, GyRelax1[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] <= G[i]);
    @constraint(smp, GyRelax2[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] <= yhat[i,par]);
    @constraint(smp, GyRelax3[i in pData.II,par in 1:length(divSet[i])], Gy[i,par] >= G[i] + yhat[i,par] - 1);

    @constraint(smp, GCons1[i in pData.II],dDω.H*G[i] - sum(H[divSet[i][par].endH]*Gy[i,par] for par in 1:length(divSet[i])) <= dDω.H - that[i]);
    @constraint(smp, GCons2[i in pData.II],sum(H[divSet[i][par].startH]*Gy[i,par] for par in 1:length(divSet[i])) - dDω.H*G[i] >= -that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i])));
    @constraint(smp, GFixed0[i in pData.II],G[i] >= sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH]));
    @constraint(smp, GFixed1[i in pData.II],G[i] <= 1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH]));

    # add the predecessors and the successors logic constraints
    @constraint(smp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(smp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(smp, tFnAnt1[i in pData.II],t[i] + G[i]*M[i] >= that[i]);
    @constraint(smp, tFnAnt2[i in pData.II],t[i] - G[i]*M[i] <= that[i]);
    @constraint(smp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(smp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(smp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(smp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(smp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(smp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(smp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(smp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(smp, Min, t[0]);
    smpStatus = solve(smp);
    vhat = getobjectivevalue(smp);
    Ghat = Dict();
    for i in pData.II
        Ghat[i] = getvalue(smp[:G][i]);
    end
    # solve the subproblem by dual formulation
    #sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, λFG1[i in pData.II, par in 1:length(divSet[i])] <= 0);
    @variable(sp, λFG2[i in pData.II, par in 1:length(divSet[i])] <= 0);
    @variable(sp, λFG3[i in pData.II, par in 1:length(divSet[i])] >= 0);
    @variable(sp, λHG1[i in pData.II] <= 0);
    @variable(sp, λHG2[i in pData.II] >= 0);
    @variable(sp, λGy1[i in pData.II] >= 0);
    @variable(sp, λGy2[i in pData.II] <= 0);
    @variable(sp, λGG[k in pData.K] <= 0);
    @variable(sp, λtG1[i in pData.II] >= 0);
    @variable(sp, λtG2[i in pData.II] >= 0);
    @variable(sp, λtG3[i in pData.II] <= 0);
    @variable(sp, λxG1[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, λxG2[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG1[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG2[i in pData.II, j in pData.Ji[i]] <= 0);
    @variable(sp, λsG3[i in pData.II, j in pData.Ji[i]] >= 0);
    @variable(sp, λbudget <= 0);
    @variable(sp, λxub[i in pData.II] <= 0);
    @variable(sp, λdur[k in pData.K] >= 0);

    tBool = Dict();
    for i in pData.II
        if i != 0
            tBool[i] = 0;
        else
            tBool[i] = 1;
        end
    end

    # t constraint
    @constraint(sp, tConstr[i in pData.II], tBool[i] - λtG1[i] - λtG2[i] - λtG3[i] +
        sum(λdur[k] for k in pData.K if k[1] == i) - sum(λdur[k] for k in pData.K if k[2] == i) >= 0);
    # x constraint
    @constraint(sp, xConstr[i in pData.II, j in pData.Ji[i]], -λxG1[i,j] - λxG2[i,j] +
        λsG2[i,j] + λsG3[i,j] - pData.b[i][j]*λbudget - λxub[i] - sum(pData.D[i]*pData.eff[i][j]*λdur[k] for k in pData.K if k[1] == i) >= 0);
    # G constraint
    @constraint(sp, gConstr[i in pData.II], sum(λFG1[i,par] + λFG3[i,par] for par in 1:length(divSet[i])) - dDω.H*(λHG1[i] - λHG2[i]) -
        λGy1[i] - λGy2[i] + sum(λGG[k] for k in pData.K if k[2] == i) - sum(λGG[k] for k in pData.K if k[1] == i) + λtG1[i]*dDω.H -
        λtG2[i]*M[i] + λtG3[i]*M[i] + sum(-λxG1[i,j] + λxG2[i,j] + λsG1[i,j] + λsG3[i,j] for j in pData.Ji[i]) +
        sum(dDω.d[i]*λdur[k] for k in pData.K if k[1] == i) >= 0);
    # Gy constraint
    @constraint(sp, FConstr[i in pData.II, par in 1:length(divSet[i])], -λFG1[i,par] - λFG2[i,par] - λFG3[i,par] +
        H[divSet[i][par].endH]*λHG1[i] - H[divSet[i][par].startH]*λHG2[i] >= 0);
    # s constraint
    @constraint(sp, sConstr[i in pData.II, j in pData.Ji[i]], -λsG1[i,j] - λsG2[i,j] - λsG3[i,j] -
        sum(dDω.d[i]*pData.eff[i][j]*λdur[k] for k in pData.K if k[1] == i) >= 0);

    @objective(sp, Max, sum(sum(yhat[i,par]*λFG2[i,par] + (yhat[i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - that[i]) + λHG2[i]*(-that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(that[i]*(λtG2[i] + λtG3[i]) + sum(xhat[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K));
    solve(sp);

    # objective function of the binary feasible solution should be the same
    @constraint(sp, binaryTight, sum(sum(yhat[i,par]*λFG2[i,par] + (yhat[i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - that[i]) + λHG2[i]*(-that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(that[i]*(λtG2[i] + λtG3[i]) + sum(xhat[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K) >= (1- 1e-6)*vhat);

    @expression(sp, hatPoint, sum(sum(yhat[i,par]*λFG2[i,par] + (yhat[i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - that[i]) + λHG2[i]*(-that[i] + sum(yhat[i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(that[i]*(λtG2[i] + λtG3[i]) + sum(xhat[i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K));

    # optimize the fractional solution's objective
    LL = length(ycoreList);
    @objective(sp, Max, sum(sum(sum(ycoreList[ll][i,par]*λFG2[i,par] + (ycoreList[ll][i,par] - 1)*λFG3[i,par] for par in 1:length(divSet[i])) for i in pData.II) +
        sum(λHG1[i]*(dDω.H - tcoreList[ll][i]) + λHG2[i]*(-tcoreList[ll][i] + sum(ycoreList[ll][i,par]*H[divSet[i][par].startH] for par in 1:length(divSet[i]))) for i in pData.II) +
        sum(λGy1[i]*(sum(ycoreList[ll][i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH])) +
            λGy2[i]*(1 - sum(ycoreList[ll][i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH])) for i in pData.II) +
        sum(tcoreList[ll][i]*(λtG2[i] + λtG3[i]) + sum(xcoreList[ll][i,j]*(λxG1[i,j] + λxG2[i,j]) for j in pData.Ji[i]) for i in pData.II) -
        sum(sum(λsG3[i,j] for j in pData.Ji[i]) for i in pData.II) + pData.B*λbudget + sum(λxub[i] for i in pData.II) +
        sum(pData.D[k[1]]*λdur[k] for k in pData.K) for ll in 1:LL));

    spStatus = solve(sp);
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    γdict = Dict();             # dual for y
    if spStatus != :Optimal
        println(ωCurr);
        hPv = NaN;
    else
        hPv = getvalue(hatPoint);

        for i in pData.II
            πdict[i] = -getvalue(sp[:λHG1][i]) - getvalue(sp[:λHG2][i]) + getvalue(sp[:λtG2][i]) + getvalue(sp[:λtG3][i]);
            for j in pData.Ji[i]
                λdict[i,j] = getvalue(sp[:λxG1][i,j]) + getvalue(sp[:λxG2][i,j]);
            end
            for par in 1:length(divSet[i])
                γdict[i,par] = H[divSet[i][par].startH]*getvalue(sp[:λHG2][i]) +
                    getvalue(sp[:λFG2][i,par]) + getvalue(sp[:λFG3][i,par]);
                if dDω.H <= H[divSet[i][par].startH]
                    γdict[i,par] += getvalue(sp[:λGy1][i]);
                elseif dDω.H >= H[divSet[i][par].endH]
                    γdict[i,par] -= getvalue(sp[:λGy2][i]);
                end
            end
        end
    end

    if returnOpt == 0
        if (spStatus == :Optimal)&(smpStatus == :Optimal)
            return πdict,λdict,γdict,hPv,Ghat;
        else
            errorHat = "Error";
            return πdict,λdict,γdict,hPv,errorHat;
        end
    else
        return πdict,λdict,γdict,hPv,sp;
    end
end

function sub_divTn(pData,dDω,ωCurr,that,xhat,yhat,divSet,H,M,returnOpt = 0)
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    # sp = Model(solver = CplexSolver(CPX_PARAM_SCRIND = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= G[i in pData.II] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    @constraint(sp, GCons1[i in pData.II],dDω.H*G[i] <= sum(H[divSet[i][par].endH]*G[i]*yhat[i,par] for par in 1:length(divSet[i])) + dDω.H - that[i]);
    @constraint(sp, GCons2[i in pData.II],-dDω.H*G[i] >= -that[i] + sum(H[divSet[i][par].startH]*(1 - G[i])*yhat[i,par] for par in 1:length(divSet[i])));
    @constraint(sp, GFixed0[i in pData.II],G[i] >= sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H <= H[divSet[i][par].startH]));
    @constraint(sp, GFixed1[i in pData.II],G[i] <= 1 - sum(yhat[i,par] for par in 1:length(divSet[i]) if dDω.H >= H[divSet[i][par].endH]));

    # add the predecessors and the successors logic constraints
    @constraint(sp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M[i] >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M[i] <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(sp, Min, t[0]);

    # obtain the dual variables for cuts
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v - λ(x - xhat) - π(t - that)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    γdict = Dict();             # dual for y
    for i in pData.II
        πdict[i] = -getdual(sp[:GCons1][i]) - getdual(sp[:GCons2][i]) +
            (getdual(sp[:tFnAnt1][i]) + getdual(sp[:tFnAnt2][i]));
        for j in pData.Ji[i]
            λdict[i,j] = (getdual(sp[:xFnAnt1][i,j]) + getdual(sp[:xFnAnt2][i,j]));
        end
        for par in 1:length(divSet[i])
            γdict[i,par] = H[divSet[i][par].endH]*getdual(sp[:GCons1][i])*getvalue(sp[:G][i]) +
                H[divSet[i][par].startH]*getdual(sp[:GCons2][i])*(1 - getvalue(sp[:G][i]));
            if dDω.H <= H[divSet[i][par].startH]
                γdict[i,par] += getdual(sp[:GFixed0][i]);
            elseif dDω.H >= H[divSet[i][par].endH]
                γdict[i,par] -= getdual(sp[:GFixed1][i]);
            end
        end
    end
    Ghat = Dict();
    for i in pData.II
        Ghat[i] = getvalue(sp[:G][i]);
    end
    if returnOpt == 0
        return πdict,λdict,γdict,vk,Ghat;
    else
        return πdict,λdict,γdict,vk,sp;
    end
end

function sub_divTDualTn(pData,dDω,ωCurr,that,xhat,yhat,divSet,H,M,tcore,xcore,ycore,returnOpt = 0)
end
