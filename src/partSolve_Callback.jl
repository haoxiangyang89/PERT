Tmax = disData[length(Ω)].H + longestPath(pData)[0];
pdData = deepcopy(pData);
for i in pData.II
    if i != 0
        pdData.D[i] = pData.D[i] + maximum([disData[ω].d[i] for ω in Ω])
    else
        pdData.D[i] = pData.D[i];
    end
end
lDict = longestPath(pdData);
for i in pData.II
    lDict[i] += disData[length(Ω)].H;
end
Tmax1 = lDict[0];

# start with an upper bound based on the deterministic solution
tdet,xdet,fdet = detBuild(pData);
ubdet = ubCalP(pData,disData,Ω,xdet,tdet,Tmax1);
brInfo = precludeRel(pData,disData,Ω,ubdet);

H = Dict();
H[0] = 0;
H[length(Ω)+1] = Tmax;
for ω in Ω
    H[ω] = disData[ω].H;
end

# initialize cutSet and divSet
cutSet = [];
divSet = Dict();
divDet = Dict();
for i in pData.II
    set1 = [ω for ω in Ω if brInfo[findfirst(pData.II,i),ω] == 1];
    setn1 = [ω for ω in Ω if brInfo[findfirst(pData.II,i),ω] == -1];

    if set1 != []
        set1t = partType(0,maximum(set1));
        if setn1 != []
            setn1t = partType(minimum(setn1),length(Ω) + 1);
            set0t = partType(maximum(set1),minimum(setn1));
            divSet[i] = [set1t,set0t,setn1t];
            divDet[i] = [1,0,-1];
        else
            set0t = partType(maximum(set1),length(Ω) + 1);
            divSet[i] = [set1t,set0t];
            divDet[i] = [1,0];
        end
    else
        if setn1 != []
            setn1t = partType(minimum(setn1),length(Ω) + 1);
            set0t = partType(0,minimum(setn1));
            divSet[i] = [set0t,setn1t];
            divDet[i] = [0,-1];
        else
            set0t = partType(0,length(Ω) + 1);
            divSet[i] = [set0t];
            divDet[i] = [0];
        end
    end
end
xbest = Dict();
tbest = Dict();
lbCost = -Inf;
lbCostList = [];
ubCostList = [ubdet];
GList = [];
cutynSel = Dict();
cutThreshold = [10];

function partBenders(cb)
    # the callback function
    that = Dict();
    xhat = Dict();
    θhat = Dict();
    yhat = Dict();
    # obtain the solution at the current node
    for i in pData.II
        that[i] = getvalue(t[i]);
        for j in pData.Ji[i]
            xhat[i,j] = getvalue(x[i,j]);
        end
        for par in 1:length(divSet[i])
            yhat[i,par] = getvalue(y[i,par]);
        end
    end
    for ω in Ω
        θhat[ω] = getvalue(θ[ω]);
    end
    # examine whether the previously generated cuts are tight
    cutSel = examineCuts_count(disData,Ω,cutSel,cutSet,that,xhat,θhat,yhat,cutThreshold[length(cutThreshold)]);

    # generate cuts
    πdict = Dict();
    λdict = Dict();
    γdict = Dict();
    vk = Dict();
    θInt = Dict();
    ubCost = minimum(ubCostList);
    ubTemp,θInt = ubCalP(pData,disData,Ω,xhat,that,Tmax1,1);
    if ubCost > ubTemp
        tbest = copy(that);
        xbest = copy(xhat);
    end
    push!(ubCostList,ubTemp);
    dataList = pmap(ω -> sub_divT(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict), Ω);
    for ω in Ω
        πdict[ω] = dataList[ω][1];
        λdict[ω] = dataList[ω][2];
        γdict[ω] = dataList[ω][3];
        vk[ω] = dataList[ω][4];
    end
    cutDual = [];
    for ω in Ω
        if vk[ω] - θhat[ω] > 1e-4*θhat[ω]
            push!(cutDual,[ω,vk[ω],πdict[ω],λdict[ω],γdict[ω]]);
            @lazyconstraint(cb, θ[ω] >= vk[ω] + sum(πdict[ω][i]*(t[i] - that[i]) for i in pData.II) +
                sum(sum(λdict[ω][i,j]*(x[i,j] - xhat[i,j]) for j in pData.Ji[i]) for i in pData.II) +
                sum(sum(γdict[ω][i,par]*(y[i,par] - yhat[i,par]) for par in 1:length(divSet[i])) for i in pData.II));
            #mp = addtxyCut(pData,ω,mp,πdict[ω],λdict[ω],γdict[ω],vk[ω],that,xhat,yhat,divSet);
        end
    end
    push!(cutSet,[[that,xhat,yhat,divSet],cutDual]);
    for l in length(cutSet[length(cutSet)])
        cutSel[length(cutSet),l] = 0;
    end
    GCurrent = [dataList[ω][5] for ω in Ω];
    push!(GList,GCurrent);
end

keepIter = true;
while keepIter
    tCurrent = Dict();
    xCurrent = Dict();
    θCurrent = Dict();
    yCurrent = Dict();

    # move the createMaster_Callback here
    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-9,FeasibilityTol = 1e-9,OutputFlag = 0));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
      y[i in pData.II, par in 1:length(divSet[i])], Bin
    end);
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(mp, tub[i in pData.II], t[i] <= sum(H[divSet[i][par].endH]*y[i,par] for par in 1:length(divSet[i])));
    @constraint(mp, tlb[i in pData.II], t[i] >= sum(H[divSet[i][par].startH]*y[i,par] for par in 1:length(divSet[i])));
    @constraint(mp, yConstr[i in pData.II], sum(y[i,par] for par in 1:length(divSet[i])) == 1);
    @constraint(mp, yLimit[i in pData.II, par in 1:length(divSet[i]); divDet[i][par] != 0], y[i,par] == 0);

    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));

    # add the cut
    # cutInfo = 2 dimensional vector, first dimention record the primal solution,
    # second dimension record the dual solution for every scenario
    for nc in 1:length(cutSet)
        for l in 1:length(cutSet[nc][2])
            if cutSel[nc,l] <= cutThreshold[length(cutThreshold)]
                ω = cutSet[nc][2][l][1];
                vk = cutSet[nc][2][l][2];
                πk = cutSet[nc][2][l][3];
                λk = cutSet[nc][2][l][4];
                γk = cutSet[nc][2][l][5];
                rhsExpr = vk;
                @constraint(mp, θ[ω] >= vk + sum(πk[i]*(mp[:t][i] - cutSet[nc][1][1][i]) +
                    sum(λk[i,j]*(mp[:x][i,j] - cutSet[nc][1][2][i,j]) for j in pData.Ji[i]) +
                    sum(γk[i,par]*(sum(mp[:y][i,parNew] for parNew in 1:length(divSet[i]) if revPar(cutSet[nc][1][4][i],divSet[i][parNew]) == par) - cutSet[nc][1][3][i,par])
                    for par in 1:length(cutSet[nc][1][4][i])) for i in pData.II));
            end
        end
    end

    # add the constraints between y
    # obtain the set of y's all predecessors
    for i in pData.II
        # for each precedence relationship
        for j in allSucc[i]
            k = (i,j);
            for par1 in 1:length(divSet[k[1]])
                for par2 in 1:length(divSet[k[2]])
                    if H[divSet[k[2]][par2].endH] < H[divSet[k[1]][par1].startH] + distanceDict[i,j]
                        @constraint(mp, y[k[1],par1] + y[k[2],par2] <= 1);
                    end
                end
            end
        end
    end

    addlazycallback(mp, partBenders);
    tic();
    solve(mp);
    tIter = toc();
    lbCost = getobjectivevalue(mp);

    # need to come up with a rule to partition: gradient descent like binary search
    # check θInt vs. θhat: why the lower bound and the upper bound do not converge quickly --->
    # use the sub problem solution G to learn the b&b
    # also need to think up a way to tightening the cuts for each partition
    ubCost = minimum(ubCostList);
    if (ubCost - lbCost)/ubCost < ϵ
        keepIter = false;
    else
        GCurrent = GList[length(GList)];
        GFrac = Dict();
        for i in pData.II
            GFraciList = [ω for ω in Ω if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)];
            if GFraciList != []
                GFrac[i] = [minimum(GFraciList),maximum(GFraciList)];
            else
                GFrac[i] = [];
            end
        end
        newPartition = [];
        for i in pData.II
            if GFrac[i] != []
                newItem = (i,Int(floor((GFrac[i][1] + GFrac[i][2])/2)));
                push!(newPartition,newItem);
            end
        end
        divSet,divDet = splitPar(divSet,divDet,newPartition);
    end
    push!(cutThreshold,cutThreshold[length(cutThreshold)] + 5);
end

# need a cut selection process within the callback