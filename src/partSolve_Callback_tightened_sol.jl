# calculate the largest possible starting time for each activity in any disrupted scenario
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

allSucc = findSuccAll(pData);
distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

H = Dict();
H[0] = 0;
# remove the duplicate ones
counter = 1;
for ω in Ω
    if !(disData[ω].H in values(H))
        H[counter] = disData[ω].H;
        counter += 1;
    end
end
H[counter] = Tmax;
HΩ = 1:(counter - 1);
HRev = Dict();
for hIter in keys(H)
    HRev[H[hIter]] = hIter;
end

# start with an upper bound based on the smaller stochastic solution
ubList,tHList,ubInc,tbest,xbest,θbest,textList,xextList = iniPart(pData,disData,Ω,sN,MM,1);
lbCost = -Inf;
lbCostList = [];
ubCost = ubInc;
ubCostList = [ubCost];

GList = [];
cutSel = Dict();
tcoreList = [];
xcoreList = [];
ycoreList = [];
errorList = [];

brInfo = precludeRelNew(pData,H,ubCost);

# initialize cutSet and divSet
cutSet = [];
divSet = Dict();
divDet = Dict();
for i in pData.II
    set1 = [h for h in HΩ if brInfo[findfirst(pData.II,i),h] == 1];
    setn1 = [h for h in HΩ if brInfo[findfirst(pData.II,i),h] == -1];

    if set1 != []
        set1t = partType(0,maximum(set1));
        if setn1 != []
            setn1t = partType(minimum(setn1),length(HΩ) + 1);
            set0t = partType(maximum(set1),minimum(setn1));
            divSet[i] = [set1t,set0t,setn1t];
            divDet[i] = [1,0,-1];
        else
            set0t = partType(maximum(set1),length(HΩ) + 1);
            divSet[i] = [set1t,set0t];
            divDet[i] = [1,0];
        end
    else
        if setn1 != []
            setn1t = partType(minimum(setn1),length(HΩ) + 1);
            set0t = partType(0,minimum(setn1));
            divSet[i] = [set0t,setn1t];
            divDet[i] = [0,-1];
        else
            set0t = partType(0,length(HΩ) + 1);
            divSet[i] = [set0t];
            divDet[i] = [0];
        end
    end
end

# pre-separate the partition
#divSet,divDet = splitPar_CI(divSet,divDet,tHList);
tcoreList = [];
xcoreList = [];
ycoreList = [];
for l in 1:length(textList)
    push!(tcoreList,textList[l]);
    push!(xcoreList,xextList[l]);
    ycoreItem = Dict();
    for i in pData.II
        for par in 1:length(divSet[i])
            if (textList[l][i] >= H[divSet[i][par].startH])&(textList[l][i] < H[divSet[i][par].endH])
                ycoreItem[i,par] = 1;
            elseif (abs(textList[l][i] - H[length(H) - 1]) < 1e-4)&(divSet[i][par].endH == length(H) - 1)
                ycoreItem[i,par] = 1;
            else
                ycoreItem[i,par] = 0;
            end
        end
    end
    push!(ycoreList,ycoreItem);
end

function partBenders(cb)
    currentLB = MathProgBase.cbgetbestbound(cb);
    println("lazy,$(currentLB)");
    if currentLB <= minimum(ubCostList)
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
                yhat[i,par] = round(getvalue(y[i,par]));
            end
        end
        for ω in Ω
            θhat[ω] = getvalue(θ[ω]);
        end
        push!(tcoreList,that);
        push!(xcoreList,xhat);
        push!(ycoreList,yhat);

        # generate cuts
        πdict = Dict();
        λdict = Dict();
        γdict = Dict();
        vk = Dict();
        θInt = Dict();
        ubCost = minimum(ubCostList);
        ubTemp,θInt = ubCalP(pData,disData,Ω,xhat,that,Tmax1,1);
        if ubCost > ubTemp
            for i in pData.II
                tbest[i] = that[i];
                for j in pData.Ji[i]
                    xbest[i,j] = xhat[i,j];
                end
            end
        end
        push!(ubCostList,ubTemp);

        #dataList = pmap(ω -> sub_divT(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict), Ω);
        # obtain the cores
        tcore,xcore,ycore = avgCore(pData,divSet,tcoreList,xcoreList,ycoreList);
        dataList = pmap(ω -> sub_divTDualT2(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore), Ω);
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
            end
        end
        push!(cutSet,[[deepcopy(that),deepcopy(xhat),deepcopy(yhat),deepcopy(divSet)],cutDual]);
        GCurrent = [dataList[ω][5] for ω in Ω];
        push!(GList,GCurrent);
    else
        return JuMP.StopTheSolver;
    end
end

keepIter = true;
lbHist = [];
ubHist = [];
timeHist = [];
cutHist = [];
intSolHist = [];
yhistList = [];
mp = Model(solver = GurobiSolver(IntFeasTol = 1e-8, FeasibilityTol = 1e-8, Threads = noThreads));
# mp = Model(solver = GurobiSolver(Threads = noThreads));
@variables(mp, begin
  θ[Ω] >= 0
  0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
  t[i in pData.II] >= 0
  y[i in pData.II, par in 1:length(divSet[i])], Bin
end);


while keepIter
    # correct all ycoreList
    for ll in 1:length(ycoreList)
        for i in pData.II
            for par in 1:length(divSet[i])
                if (tcoreList[ll][i] >= H[divSet[i][par].startH])&(tcoreList[ll][i] < H[divSet[i][par].endH])
                    ycoreList[ll][i,par] = 1;
                else
                    ycoreList[ll][i,par] = 0;
                end
            end
        end
    end

    # move the createMaster_Callback here
    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-8, FeasibilityTol = 1e-8, Threads = noThreads));
    # mp = Model(solver = GurobiSolver(Threads = noThreads));
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

    #@objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));
    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω) + r*sum(t[i] for i in pData.II));

    tCurrent = Dict();
    xCurrent = Dict();
    θCurrent = Dict();
    yCurrent = Dict();

    # add the cut
    # cutInfo = 2 dimensional vector, first dimention record the primal solution,
    # second dimension record the dual solution for every scenario
    for nc in 1:length(cutSet)
        for l in 1:length(cutSet[nc][2])
            ω = cutSet[nc][2][l][1];
            vk = cutSet[nc][2][l][2];
            πk = cutSet[nc][2][l][3];
            λk = cutSet[nc][2][l][4];
            γk = cutSet[nc][2][l][5];
            @constraint(mp, mp[:θ][ω] >= vk + sum(πk[i]*(mp[:t][i] - cutSet[nc][1][1][i]) +
                sum(λk[i,j]*(mp[:x][i,j] - cutSet[nc][1][2][i,j]) for j in pData.Ji[i]) +
                sum(γk[i,par]*(sum(mp[:y][i,parNew] for parNew in 1:length(divSet[i]) if revPar(cutSet[nc][1][4][i],divSet[i][parNew]) == par) - cutSet[nc][1][3][i,par])
                for par in 1:length(cutSet[nc][1][4][i])) for i in pData.II));
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

    # find the current best upper bound from coreList
    # lbmin = ubInc;
    # coreMinind = -1;
    # θbestR = Dict();
    # for l in 1:length(tcoreList)
    #     θcore = pmap(ω -> sub_divT(pData,disData[ω],ω,tcoreList[l],xcoreList[l],ycoreList[l],divSet,H,lDict),Ω);
    #     lbcore = pData.p0*(tcoreList[l][0]) + sum(θcore[ω]*disData[ω].prDis for ω in Ω);
    #     if lbcore < lbmin
    #         lbmin = lbcore;
    #         coreMinind = l;
    #         for ω in 1:length(Ω)
    #             θbestR[Ω[ω]] = θcore[ω];
    #         end
    #     end
    # end
    ybest = Dict();
    for i in pData.II
        setvalue(t[i], tbest[i]);
        for j in pData.Ji[i]
            setvalue(x[i,j],xbest[i,j]);
        end
        for par in 1:length(divSet[i])
            if (tbest[i] >= H[divSet[i][par].startH])&(tbest[i] < H[divSet[i][par].endH])
                ybest[i,par] = 1;
            elseif (abs(tbest[i] - H[length(H) - 1]) < 1e-4)&(divSet[i][par].endH == length(H) - 1)
                ybest[i,par] = 1;
            else
                ybest[i,par] = 0;
            end
            setvalue(y[i,par],ybest[i,par]);
        end
    end
    θbest = pmap(ω -> sub_divT(pData,disData[ω],ω,tbest,xbest,ybest,divSet,H,lDict),Ω);
    for ω in Ω
        setvalue(θ[ω],θbest[ω]);
    end

    addlazycallback(mp, partBenders);
    #global setSolBool = true;
    #addheuristiccallback(mp, setSolcb);
    tic();
    solve(mp);
    tIter = toc();
    push!(timeHist,tIter);

    lbCost = getobjectivevalue(mp);
    push!(lbHist,lbCost);

    # only select currently tight cuts
    for i in pData.II
        tCurrent[i] = getvalue(mp[:t][i]);
        for j in pData.Ji[i]
            xCurrent[i,j] = getvalue(mp[:x][i,j]);
        end
        for par in 1:length(divSet[i])
            yCurrent[i,par] = getvalue(mp[:y][i,par]);
        end
    end
    for ω in Ω
        θCurrent[ω] = getvalue(mp[:θ][ω]);
    end

    # need to come up with a rule to partition: gradient descent like binary search
    # check θInt vs. θhat: why the lower bound and the upper bound do not converge quickly --->
    # use the sub problem solution G to learn the b&b
    # also need to think up a way to tightening the cuts for each partition
    ubCost = minimum(ubCostList);
    push!(ubHist,ubCost);
    push!(intSolHist,length(ubCostList));
    if (ubCost - lbCost)/ubCost < ϵ
        keepIter = false;
    else
        # cutSel = examineCuts_count_2(disData,Ω,cutSet,divSet,tCurrent,xCurrent,θCurrent,yCurrent);
        # cutSetNew = selectCuts2(cutSet,cutSel);
        GCurrent = GList[length(GList)];
        GFrac = Dict();
        for i in pData.II
            GFraciList = [disData[ω].H for ω in Ω if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)];
            if GFraciList != []
                GFrac[i] = [HRev[minimum(GFraciList)],HRev[maximum(GFraciList)]];
            else
                GFrac[i] = [];
            end
        end
        newPartition = [];
        for i in pData.II
            if GFrac[i] != []
                newItem = (i,GFrac[i][1],GFrac[i][2]);
                push!(newPartition,newItem);
            end
        end
        #divSet,divDet = splitPar(divSet,divDet,newPartition);
        divSet,divDet = splitPar3(divSet,divDet,newPartition);
        push!(cutHist,sum(length(cutSet[l][2]) for l in 1:length(cutSet)));
        # cutSet = deepcopy(cutSetNew);
    end
end

# need a cut selection process within the callback
