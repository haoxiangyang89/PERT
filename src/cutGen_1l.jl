# This is the script with all cut generation
function createMaster(pData,disData,Ω)
    # mp = Model(solver = GurobiSolver(OutputFlag = 0));
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
    end);
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @objective(mp,Min,pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    return mp;
end

function bGenbuild(pData,dDω,xhat,that,brInfo)
    M = sum(max(pData.D[i],pData.D[i]+dDω.d[i]) for i in pData.II if i != 0);

    # sp = Model(solver = GurobiSolver(OutputFlag = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= F[i in pData.II] <= 1);
    @variable(sp, 0 <= G[i in pData.II] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0],dDω.H - F[i]*M <= that[i]);
    @constraint(sp, GCons[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0],dDω.H - 1e-6 + G[i]*M >= that[i]);
    @constraint(sp, FFixed[i in pData.II; brInfo[findin(pData.II,i)[1]] == -1],F[i] == 1);
    @constraint(sp, GFixed[i in pData.II; brInfo[findin(pData.II,i)[1]] == 1],G[i] == 1);
    @constraint(sp, FGCons[i in pData.II],F[i] + G[i] == 1);

    # add the predecessors and the successors logic constraints
    @constraint(sp, FPredecessors[i in pData.II, k in pData.Pre[i]], F[i] <= F[k]);
    @constraint(sp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + (1 - F[i])*M >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - (1 - F[i])*M <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + (1 - F[i]) >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - (1 - F[i]) <= xhat[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j] + 1 - G[i]);
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

function bGenbuild_New(pData,dDω,xhat,that,brInfoω,M = Dict())
    # imbed the logic in the program
    if M == Dict()
        MM = iSolve_NC(pData,dDω,0,brInfoω);
        for i in pData.II
            M[i] = MM;
        end
    end

    # sp = Model(solver = GurobiSolver(OutputFlag = 0));
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= F[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0] <= 1);
    @variable(sp, 0 <= G[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]; brInfoω[findin(pData.II,i)[1]] == 0] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0],dDω.H - F[i]*M[i] <= that[i]);
    @constraint(sp, GCons[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0],dDω.H + G[i]*M[i] >= that[i]);
    @constraint(sp, FGCons[i in pData.II; brInfoω[findin(pData.II,i)[1]] == 0],F[i] + G[i] == 1);

    # add the predecessors and the successors logic constraints
    @constraint(sp, FPredecessors[i in pData.II, k in pData.Pre[i]; (brInfoω[findin(pData.II,i)[1]] == 0)&(brInfoω[findin(pData.II,k)[1]] == 0)], F[i] <= F[k]);
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

function subInt(pData,dDω,xhat,that)
    # solve the MIP recourse problem
    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);

    # add the basic sub problem constraints
    @constraint(sp, fixT[i in pData.II; that[i] < dDω.H], t[i] == that[i]);
    @constraint(sp, boundT[i in pData.II; that[i] >= dDω.H], t[i] >= dDω.H);
    @constraint(sp, fixX[i in pData.II,j in pData.Ji[i]; that[i] < dDω.H], x[i,j] == xhat[i,j]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K; that[k[1]] < dDω.H], t[k[2]] - t[k[1]] >= pData.D[k[1]]*
        (1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(sp, durationConstr2[k in pData.K; that[k[1]] >= dDω.H], t[k[2]] - t[k[1]] >= (pData.D[k[1]]+dDω.d[k[1]])*
        (1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));

    @objective(sp, Min, t[0]);
    solve(sp);
    return getobjectivevalue(sp);
end

function ubCal(pData,disData,Ω,xhat,that)
    # calculate the upper bound of the problem given the master solution
    ubCost = that[0]*pData.p0;
    for ω in Ω
        cω = subInt(pData,disData[ω],xhat,that);
        ubCost += disData[ω].prDis*cω;
    end
    return ubCost;
end

function ubCalP(pData,disData,Ω,xhat,that)
    # parallel version of calculating the upper bound
    ubCost = that[0]*pData.p0;
    cωList = pmap(ω -> subInt(pData,disData[ω],xhat,that), Ω);
    ubCost += sum(cωList[i]*disData[Ω[i]].prDis for i in 1:length(ω));
    return ubCost;
end

# build the LP to obtain the earliest possible starting time of an activity
function iSolve(pData,iTarget)
    # mp = Model(solver = GurobiSolver(OutputFlag = 0));
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(mp, t[i in pData.II] >= 0);
    @variable(mp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);

    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 -
        sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, budgetConstr, sum(sum(x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    @objective(mp, Min, t[iTarget]);
    solve(mp);

    return getobjectivevalue(mp);
end

# build the LP to obtain the big M
function iSolve_NC(pData,dDω,iTarget,brInfoω)
    # mp = Model(solver = GurobiSolver(OutputFlag = 0));
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(mp, t[i in pData.II] >= 0);
    Dω = Dict();

    for i in keys(dDω.d)
        if brInfoω[findin(pData.II,i)[1]] == 0
            Dω[i] = max(0,dDω.d[i])+pData.D[i];
        elseif brInfoω[findin(pData.II,i)[1]] == -1
            Dω[i] = pData.D[i];
        else
            Dω[i] = pData.D[i] + dDω.d[i];
        end
    end
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= Dω[k[1]]);

    @objective(mp, Min, t[iTarget]);
    solve(mp);

    return getobjectivevalue(mp);
end

function bGenbuild_Dual(pData,dDω,xhat,that,brInfoω,M = Dict())
    # imbed the logic in the program
    if M == Dict()
        MM = iSolve_NC(pData,dDω,0,brInfoω);
        for i in pData.II
            M[i] = MM;
        end
    end

    dp = Model(solver = GurobiSolver(OutputFlag = 0));

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
