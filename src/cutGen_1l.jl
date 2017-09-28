# This is the script with all cut generation
function createMaster(pData,disData,Ω)
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

function bGenbuild_New(pData,dDω,xhat,that,brInfo)
    # imbed the logic in the program
    M = sum(max(pData.D[i],pData.D[i]+dDω.d[i]) for i in pData.II if i != 0);

    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= F[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0] <= 1);
    @variable(sp, 0 <= G[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]; brInfo[findin(pData.II,i)[1]] == 0] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0],dDω.H - F[i]*M <= that[i]);
    @constraint(sp, GCons[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0],dDω.H - 1e-6 + G[i]*M >= that[i]);
    @constraint(sp, FGCons[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0],F[i] + G[i] == 1);

    # add the basic sub problem constraints for the undecided activities
    @constraint(sp, tGbound1[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0],t[i] + (1 - F[i])*M >= that[i]);
    @constraint(sp, tFnAnt2[i in pData.II; brInfo[findin(pData.II,i)[1]] == 0],t[i] - (1 - F[i])*M <= that[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1]] == 0],x[i,j] + (1 - F[i])*M >= xhat[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1]] == 0],x[i,j] - (1 - F[i])*M <= xhat[i,j]);

    # add the logic constraints for the activities that have already started
    @constraint(sp, tFnAnt3[i in pData.II; brInfo[findin(pData.II,i)[1]] == -1], t[i] == that[i]);
    @constraint(sp, xFnAnt3[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1]] == -1],x[i,j] == xhat[i,j]);
    @constraint(sp, tGbound2[i in pData.II; brInfo[findin(pData.II,i)[1]] == 1],t[i] >= dDω.H);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1]] == 0], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1]] == 0], s[i,j] <= x[i,j] + 1 - G[i]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]; brInfo[findin(pData.II,i)[1]] == 0], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K; brInfo[findin(pData.II,k[1])[1]] == 0], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));
    @constraint(sp, durationConstr2[k in pData.K; brInfo[findin(pData.II,k[1])[1]] == 1], t[k[2]] - t[k[1]] >= (pData.D[k[1]] + dDω.d[k[1]])*(1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(sp, durationConstr3[k in pData.K; brInfo[findin(pData.II,k[1])[1]] == -1], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));

    @objective(sp, Min, t[0]);

    # obtain the dual variables for cuts
    solve(sp);
    vk = getobjectivevalue(sp);
    # the cut generated is θ >= v - λ(x - xhat) - π(t - that)
    λdict = Dict();             # dual for x
    πdict = Dict();             # dual for t
    for i in pData.II
        if brInfo[findin(pData.II,i)[1]] == 0
            πdict[i] = -(getdual(sp[:FCons])[i] + getdual(sp[:GCons])[i] + getdual(sp[:tFnAnt1])[i] + getdual(sp[:tFnAnt2])[i]);
            for j in pData.Ji[i]
                λdict[i,j] = -(getdual(sp[:xFnAnt1])[i,j] + getdual(sp[:xFnAnt2])[i,j]);
            end
        elseif brInfo[findin(pData.II,i)[1]] == -1
            πdict[i] = -(getdual(sp[:tFnAnt3])[i]);
            for j in pData.Ji[i]
                λdict[i,j] = -(getdual(sp[:xFnAnt3])[i,j]);
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
    M = sum(max(pData.D[i],pData.D[i]+dDω.d[i]) for i in pData.II if i != 0);

    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);

    # add the basic sub problem constraints
    @constraint(sp, fixT[i in pData.II; that[i] <= dDω.H - 1e-6], t[i] == that[i]);
    @constraint(sp, boundT[i in pData.II; that[i] > dDω.H - 1e-6], t[i] >= dDω.H);
    @constraint(sp, fixX[i in pData.II,j in pData.Ji[i]; that[i] < dDω.H - 1e-6], x[i,j] == xhat[i,j]);

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
function iSolve(pData,disData,iTarget)
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
