@everywhere include("genDisjunctive.jl");

# solve the linear relaxation with disjunctive cuts
function solveLR(pData,dDω,cutSetω,tm,xm,M,returnDual = 0)
    sp = Model(solver = GurobiSolver(OutputFlag = 0));

    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= G[i in pData.II] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II],dDω.H - (1 - G[i])*M <= tm[i]);
    @constraint(sp, GCons[i in pData.II],dDω.H + G[i]*M >= tm[i]);

    # add the predecessors and the successors logic constraints
    @constraint(sp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M >= tm[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M <= tm[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xm[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xm[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    cuts = Dict();
    for nc in 1:length(cutSetω)
        π0c = cutSetω[nc][1];
        λ0c = cutSetω[nc][2];
        πc = cutSetω[nc][3];
        λc = cutSetω[nc][4];
        γc = cutSetω[nc][5];
        νc = cutSetω[nc][6];
        vc = cutSetω[nc][7];
        # add the disjunctive cuts here
        cuts[nc] = @constraint(sp, vc - sum(πc[i]*t[i] + γc[i]*G[i] + sum(λc[i,j]*x[i,j] + νc[i,j]*s[i,j] for j in pData.Ji[i]) for i in pData.II)
            <= sum(π0c[i]*tm[i] + sum(λ0c[i,j]*xm[i,j] for j in pData.Ji[i]) for i in pData.II));
    end

    @objective(sp, Min, t[0]);

    solve(sp);
    v = getobjectivevalue(sp);

    if returnDual == 0
        # return the primal solution
        tDict = Dict();
        xDict = Dict();
        gDict = Dict();
        sDict = Dict();
        for i in pData.II
            tDict[i] = getvalue(sp[:t][i]);
            gDict[i] = getvalue(sp[:G][i]);
            for j in pData.Ji[i]
                xDict[i,j] = getvalue(sp[:x][i,j]);
                sDict[i,j] = getvalue(sp[:s][i,j]);
            end
        end
        return tDict,xDict,gDict,sDict,v,sp;
    else
        # return the dual solution
        πDict = Dict();
        λDict = Dict();
        for i in pData.II
            πDict[i] = getdual(FCons[i]) + getdual(GCons[i]) + getdual(tFnAnt1[i]) + getdual(tFnAnt2[i]) +
                sum(getdual(cuts[nc])*cutSetω[nc][1][i] for nc in 1:length(cutSetω));
            for j in pData.Ji[i]
                λDict[i,j] = getdual(xFnAnt1[i,j]) + getdual(xFnAnt2[i,j]) +
                    sum(getdual(cuts[nc])*cutSetω[nc][2][i,j] for nc in 1:length(cutSetω));
            end
        end
        return πDict,λDict,v,sp;
    end
end

# solve the linear relaxation with binary restrictions from the tree
function solveLR01(pData,dDω,cutSetω,tm,xm,zeroSet,oneSet,M)
    sp = Model(solver = GurobiSolver(OutputFlag = 0));

    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    # relax the logic binary variables
    @variable(sp, 0 <= G[i in pData.II] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);

    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II],dDω.H - (1 - G[i])*M <= tm[i]);
    @constraint(sp, GCons[i in pData.II],dDω.H + G[i]*M >= tm[i]);

    # add the predecessors and the successors logic constraints
    @constraint(sp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the basic sub problem constraints
    @constraint(sp, tGbound[i in pData.II],t[i] >= dDω.H*G[i]);
    @constraint(sp, tFnAnt1[i in pData.II],t[i] + G[i]*M >= tm[i]);
    @constraint(sp, tFnAnt2[i in pData.II],t[i] - G[i]*M <= tm[i]);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],x[i,j] + G[i] >= xm[i,j]);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],x[i,j] - G[i] <= xm[i,j]);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    # G binary constraints
    @constraint(sp, GBin0[i in pData.II; i in zeroSet], G[i] == 0);
    @constraint(sp, GBin1[i in pData.II; i in oneSet], G[i] == 1);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    cuts = Dict();
    for nc in 1:length(cutSetω)
        π0c = cutSetω[nc][1];
        λ0c = cutSetω[nc][2];
        πc = cutSetω[nc][3];
        λc = cutSetω[nc][4];
        γc = cutSetω[nc][5];
        νc = cutSetω[nc][6];
        vc = cutSetω[nc][7];
        # add the disjunctive cuts here
        cuts[nc] = @constraint(sp, vc - sum(πc[i]*t[i] + γc[i]*G[i] + sum(λc[i,j]*x[i,j] + νc[i,j]*s[i,j] for j in pData.Ji[i]) for i in pData.II)
            <= sum(π0c[i]*tm[i] + sum(λ0c[i,j]*xm[i,j] for j in pData.Ji[i]) for i in pData.II));
    end

    @objective(sp, Min, t[0]);

    spst = solve(sp);

    tDict = Dict();
    xDict = Dict();
    gDict = Dict();
    sDict = Dict();
    if spst == :Optimal
        v = getobjectivevalue(sp);
        for i in pData.II
            tDict[i] = getvalue(sp[:t][i]);
            gDict[i] = getvalue(sp[:G][i]);
            for j in pData.Ji[i]
                xDict[i,j] = getvalue(sp[:x][i,j]);
                sDict[i,j] = getvalue(sp[:s][i,j]);
            end
        end
    else
        v = Inf;
    end
    return tDict,xDict,gDict,sDict,v,spst;
end

function maxGfrac(pData,G)
    maxFrac = 1;
    fracInd = -1;
    for i in pData.II
        fracI = max(G[i], 1 - G[i]);
        if fracI < maxFrac
            maxFrac = fracI;
            fracInd = i
        end
    end
    return maxFrac,fracInd;
end

function middleGfrac(pData,G)
    GList = []
    for i in pData.II
        if (G[i] > 1e-5)&(G[i] < 1 - 1e-5)
            push!(GList,(i,G[i]));
        end
    end
    GList = sort(GList, by=x->x[2]);
    if isodd(length(GList))
        loc = Int64((length(GList) + 1)/2);
    else
        loc = Int64(length(GList)/2);
    end
    middleFrac = GList[loc][2];
    fracInd = GList[loc][1];

    return middleFrac,fracInd;
end

# create the B&B tree for the subproblem
function BBprocess(pData,dDω,cutSetω,tm,xm,nTree,M)
    tree = [];
    ts,xs,gs,ss,vs,sp = solveLR(pData,dDω,cutSetω,tm,xm,M,0);
    nID = 1;
    nodeLP = treeNode(nID,vs,[ts,xs,gs],-1,[],[],[]);
    push!(tree,nodeLP);
    activeNode = [nodeLP];
    UB = Inf;
    # while it is within the size
    while (length(tree) <= nTree)&(activeNode != [])
        # breadth first search
        currentNode = activeNode[1];
        shift!(activeNode);

        gBin = true;
        for i in pData.II
            if (abs(currentNode.sol[3][i] - 1) > 1e-5)&(abs(currentNode.sol[3][i]) > 1e-5)
                gBin = false;
            end
        end
        if gBin
            # if it is binary, update the upper bound
            if currentNode.objV < UB
                UB = currentNode.objV;
            end
        else
            # if it is not binary, split the node into two
            if currentNode.objV <= UB
                # find the entry to split
                # maxFrac,fracInd = maxGfrac(pData,currentNode.sol[3]);
                middleFrac,fracInd = middleGfrac(pData,currentNode.sol[3]);

                # split to two nodes
                node1zeroSet = copy(currentNode.zeroSet);
                node1oneSet = copy(currentNode.oneSet);
                push!(node1zeroSet,fracInd);
                ts1,xs1,gs1,ss1,vs1,spst1 = solveLR01(pData,dDω,cutSetω,tm,xm,node1zeroSet,node1oneSet,M);
                if spst1 == :Optimal
                    nID += 1;
                    node1 = treeNode(nID,vs1,[ts1,xs1,gs1,ss1],currentNode.nodeID,[],node1zeroSet,node1oneSet);
                    push!(currentNode.childSet,nID);
                    push!(tree,node1);
                    push!(activeNode,node1);
                end

                node2zeroSet = copy(currentNode.zeroSet);
                node2oneSet = copy(currentNode.oneSet);
                push!(node2oneSet,fracInd);
                ts2,xs2,gs2,ss2,vs2,spst2 = solveLR01(pData,dDω,cutSetω,tm,xm,node2zeroSet,node2oneSet,M);
                if spst2 == :Optimal
                    nID += 1;
                    node2 = treeNode(nID,vs2,[ts2,xs2,gs2,ss2],currentNode.nodeID,[],node2zeroSet,node2oneSet);
                    push!(currentNode.childSet,nID);
                    push!(tree,node2);
                    push!(activeNode,node2);
                end
            end
        end
    end
    leafNodes = [];
    for item in tree
        if item.childSet == []
            push!(leafNodes,item);
        end
    end
    return leafNodes;
end

# use the sub-B&B tree to update cuts
function updateCut(pData,dDω,cutSetω,leafNodes,tm,xm,M,Mt)
    # construct the disjunctive cuts from the B&B tree
    inSet = false;
    while !inSet
        # obtain the current sub solution with generated cuts
        ts,xs,gs,ss,vs,sps = solveLR(pData,dDω,cutSetω,tm,xm,M);
        # while the current solution is not within disjunctive set
        vv,viov,πv,λv,γv,νv,π0v,λ0v = genDisjunctive(pData,dDω,cutSetω,leafNodes,tm,xm,ts,xs,gs,ss,M,Mt);
        if viov < 1e-4
            inSet = true;
        else
            push!(cutSetω,(π0v,λ0v,πv,λv,γv,νv,vv));
        end
    end
    return cutSetω;
end

function testCutM(pData,vc,πc,λc,tm,xm,ttest,xtest)
    rhsv = vc + sum(πc[i]*(ttest[i] - tm[i]) for i in pData.II);
    for i in pData.II
        for j in pData.Ji[i]
            rhsv += λc[i,j]*(xtest[i,j] - xm[i,j]);
        end
    end
    return rhsv;
end

function testCutS(pData,cutSetω,tm,xm,ts,xs,gs,ss)
    for nc in 1:length(cutSetω)
       rhsv = sum(cutSetω[nc][1][i]*tm[i] + cutSetω[nc][3][i]*ts[i] + cutSetω[nc][5][i]*gs[i] for i in pData.II);
       for i in pData.II
           for j in pData.Ji[i]
               rhsv += cutSetω[nc][2][i,j]*xm[i,j] + cutSetω[nc][4][i,j]*xs[i,j] + cutSetω[nc][6][i,j]*ss[i,j];
           end
       end
       if cutSetω[nc][7] > rhsv
           println(nc," ",cutSetω[nc][7]," ",rhsv);
       end
   end
end

# convexification of the subproblem
function convexify(pData,disData,Ω,Tmax,Tmax1,nTree,ϵ)
    # while it has not reached the optimum
    stopBool = false;
    masterCuts = [];
    mp = createMaster(pData,disData,Ω,Tmax);
    LB = -Inf;
    UB = Inf;
    xbest = Dict();
    tbest = Dict();
    cutSet = Dict();
    # initialize with an empty cut set
    for ω in Ω
        cutSet[ω] = [];
    end

    while !(stopBool)
        # solve the master problem and obtain the master solution
        solve(mp);
        tm = Dict();
        xm = Dict();
        θm = Dict();
        innerEnd = Dict();
        for i in pData.II
            tm[i] = getvalue(mp[:t][i]);
            for j in pData.Ji[i]
                xm[i,j] = getvalue(mp[:x][i,j]);
            end
        end
        for ω in Ω
            θm[ω] = getvalue(mp[:θ][ω]);
            innerEnd[ω] = false;
        end

        # update the lower bound and the upper bound
        if getobjectivevalue(mp) > LB
            LB = getobjectivevalue(mp);
        end
        ubTemp,θInt = ubCalP(pData,disData,Ω,xm,tm,Tmax1,1);
        if ubTemp < UB
            UB = ubTemp;
        end

        # check if the stopping criterion is met
        if (UB - LB)/UB > ϵ
            # for each scenario, generate disjunctive cuts using BB-D algorithm
            for ω in Ω
                dDω = disData[ω];
                cutSetω = cutSet[ω];
                # if the solution is not binary, obtain a B&B tree and add disjunctive cuts
                leafω = BBprocess(pData,dDω,cutSetω,tm,xm,nTree,Tmax1);
                cutSet[ω] = updateCut(pData,dDω,cutSetω,leafω,tm,xm,Tmax1,Tmax);
                # generate master cuts
                πlr,λlr,vlr,slr = solveLR(pData,dDω,cutSetω,tm,xm,Tmax1,1);
                if θm[ω] < vlr - 1e-4
                    mp = addtxCut(pData,ω,mp,πlr,λlr,vlr,tm,xm);
                end
            end
        else
            stopBool = true;
            # record the best solution
            tbest = copy(tm);
            xbest = copy(xm);
        end
    end
    return tbest,xbest,LB,UB;
end
