function recoverInfo(II,D,eff,b,K,p0,B)
    IIList = [i for i in II];
    Ddict = Dict();
    for i in 1:length(IIList)
        Ddict[IIList[i]] = D[i];
    end
    KList = [(K[k,1],K[k,2]) for k in 1:size(K)[1]];
    Succ = Dict();
    Pre = Dict();
    for i in IIList
        Succ[i] = [];
        Pre[i] = [];
    end
    for k in 1:size(K)[1]
        push!(Succ[K[k,1]],K[k,2]);
        push!(Pre[K[k,2]],K[k,1]);
    end
    bdict = Dict();
    effdict = Dict();
    Ji = Dict();
    Jmax = size(eff)[2];
    for i in 1:length(IIList)
        Ji[IIList[i]] = [j for j in 1:Jmax if eff[i,j] != 0];
        bdict[IIList[i]] = [b[i,j] for j in 1:Jmax if b[i,j] != 0];
        effdict[IIList[i]] = [eff[i,j] for j in 1:Jmax if eff[i,j] != 0];
    end
    pd = pInfo(IIList, Ji, Ddict, bdict, effdict, B, p0, KList, Pre, Succ)
    return pd;
end

function recoverdInfo(II,HOriShare,disdShare,disPrShare)
    Ω = [];
    disData = Dict();
    for ω in 1:length(HOriShare)
        push!(Ω,ω);
        d = Dict();
        for i in 1:length(II)
            d[II[i]] = disdShare[i,ω];
        end
        disData[ω] = disInfo(HOriShare[ω],d,disPrShare[ω]);
    end
    H = unique([HOriShare[ω] for ω in Ω]);
    unshift!(H,0);
    push!(H,Tmax);
    return disData,H,Ω;
end

function recoverdis(distanceShare)
    allSucc = Dict();
    distanceDict = Dict();
    for i in 1:size(distanceShare)[1]
        if distanceShare[i,1] in keys(allSucc)
            push!(allSucc[distanceShare[i,1]],distanceShare[i,2]);
        else
            allSucc[distanceShare[i,1]] = [distanceShare[i,2]];
        end
        distanceDict[distanceShare[i,1],distanceShare[i,2]] = distanceShare[i,3];
    end
    return allSucc,distanceDict;
end

function recoverlDict(II,lDictShare)
    lDict = Dict();
    for i in 1:length(II)
        lDict[II[i]] = lDictShare[i];
    end
    return lDict;
end

function recoverDiv(divInfoShare)
    divSet = Dict();
    divDet = Dict();
    for di in 1:length(divInfoShare)
        if !(divInfoShare[di,1] in keys(divSet))
            divSet[i] = [partType(divInfoShare[di,3],divInfoShare[di,4])];
            divDet[i] = [divInfoShare[di,5]];
        else
            push!(divSet[i],partType(divInfoShare[di,3],divInfoShare[di,4]));
            push!(divDet[i],divInfoShare[di,5]);
        end
    end
    return divSet,divDet;
end

function recoverCoreList(IIShare,IJPair,textShare,xextShare,ubextShare)
    textList = [];
    xextList = [];
    ubextList = [];
    for it in 1:size(textShare)[2]
        tDict = Dict();
        for i in 1:length(IIShare)
            tDict[IIShare[i]] = textShare[i,it];
        end
        push!(textList,tDict);
    end
    for ijt in 1:size(xextShare)[2]
        xDict = Dict();
        for ij in 1:length(IJPair)
            xDict[IJPair[ij]] = xextShare[ij,ijt];
        end
        push!(xextList,xDict);
    end
    for iu in 1:length(ubextShare)
        push!(ubextList,ubextShare[iu]);
    end
    return textList,xextList,ubextList;
end

function subPara1(pData,disData,Ω,tbest,xbest,ybest,divSet,H,lDict,wp = CachingPool(workers()))
    θList = pmap(wp,ω -> sub_divT(pData,disData[ω],ω,tbest,xbest,ybest,divSet,H,lDict),Ω);
    return θList;
endD

function subPara(pData,disData,Ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore,wp = CachingPool(workers()))
    dataList = pmap(wp,ω -> sub_divTDualT2(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore), Ω);
    return dataList;
end

function testFeas(pData,H,divSet,divDet,tcoreList,ubcoreList)
    # test the feasibility of each of the past solution
    tcoreInd = -1;
    ubcoreMin = Inf;
    for tSoli in 1:length(tcoreList)
        tSol = tcoreList[tSoli];
        feasBool = true;
        for i in pData.II
            ibSet = [l for l in 1:length(divSet[i]) if divDet[i][l] == 0];
            iub = divSet[i][maximum(ibSet)].endH;
            ilb = divSet[i][maximum(ibSet)].startH;
            if (tSol[i] < ilb)|(tSol[i] >= iub)
                feasBool = false;
            end
        end
        if feasBool
            if ubcoreList[tSoli] < ubcoreMin
                tcoreInd = tSoli;
                ubcoreMin = ubcoreList[tSoli];
            end
        end
    end
    return tcoreInd;
end

function solveMP_para_Share(data)
    # input: selectNode,tcore,xcore,weigthCore,ubCost,tbest,xbest,noTh,wpList
    divSet,divDet = recoverDiv(cutList[data[1]][8];
    ancestorList = treeList[data[1]][2];
    weightCore = data[4];
    ubCost = data[5];
    tbest = data[6];
    xbest = data[7];
    noTh = data[8];
    wp = CachingPool(data[9]);

    Tmax1 =lDict[0];
    GList = [];
    tcoreNew = [];
    xcoreNew = [];
    ubcoreNew = [];

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
            push!(tcoreNew,that);
            push!(xcoreNew,xhat);

            # generate cuts
            πdict = Dict();
            λdict = Dict();
            γdict = Dict();
            vk = Dict();
            θInt = Dict();
            ubCost = minimum(ubCostList);
            ubTemp,θInt = ubCalP(pData,disData,Ω,xhat,that,Tmax1,1,wp);
            push!(ubcoreList,ubTemp);
            push!(ubcoreNew,ubTemp);

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
            # here is the issue, pack it in a function prevent separating it
            dataList = subPara(pData,disData,Ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore,wp);
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
            push!(cutSet,[[that,xhat,yhat,divSet],cutDual]);
            GCurrent = [dataList[ω][5] for ω in Ω];
            push!(GList,GCurrent);
        else
            return JuMP.StopTheSolver;
        end
    end

    # correct all ycoreList
    ubCostList = [ubCost];
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
    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-8, FeasibilityTol = 1e-8, Threads = noTh));
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

    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));
    # @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω) + sum((t[i] - tbest[i])^2 for i in pData.II));

    # find a feasible solution and plug in the best feasible solution
    tcoreInd = testFeas(pData,H,divSet,divDet,tcoreList,ubcoreList);
    if tcoreInd != -1
        yfeas = Dict();
        tfeas = tcoreList[tcoreInd];
        xfeas = xcoreList[xcoreInd];
        for i in pData.II
            setvalue(t[i], tfeas[i]);
            for j in pData.Ji[i]
                setvalue(x[i,j],xfeas[i,j]);
            end
            for par in 1:length(divSet[i])
                if (tfeas[i] >= H[divSet[i][par].startH])&(tfeas[i] < H[divSet[i][par].endH])
                    yfeas[i,par] = 1;
                elseif (abs(tfeas[i] - H[length(H) - 1]) < 1e-4)&(divSet[i][par].endH == length(H) - 1)
                    yfeas[i,par] = 1;
                else
                    yfeas[i,par] = 0;
                end
                setvalue(y[i,par],yfeas[i,par]);
            end
        end
        θfeas = subPara1(pData,disData,Ω,tfeas,xfeas,yfeas,divSet,H,lDict,wp);
        for ω in Ω
            setvalue(θ[ω],θfeas[ω]);
        end
    end

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
            @constraint(mp, θ[ω] >= vk + sum(πk[i]*(mp[:t][i] - cutSet[nc][1][1][i]) +
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
                    if H[divSet[k[2]][par2].endH] < H[divSet[k[1]][par1].startH] + distanceDict[k[1],k[2]]
                        @constraint(mp, y[k[1],par1] + y[k[2],par2] <= 1);
                    end
                end
            end
        end
    end
    addlazycallback(mp, partBenders);
    mpStatus = solve(mp);
    mpObj = getobjectivevalue(mp);

    return mpStatus,mpObj,GList,tbest,xbest,minimum(ubCostList),cutSet,tcoreNew,xcoreNew,ycoreNew,ubcoreNew;
end

function runPara_Share(treeList,cutList,tcoreList,xcoreList,ubCost,tbest,xbest,batchNo)
    npList = [ib for ib in 2:batchNo+1];
    wpList = [ib for ib in workers() if !(ib in npList)];
    global keepIter = true;
    global noTh = div(noThreads,batchNo);

    @sync begin
        for ip in 1:length(npList)
            p = npList[ip];
            @async begin
                while true
                    # if all nodes are processed and no nodes are being processed, exit
                    boolFinished = true;
                    for l in 1:length(treeList)
                        if treeList[l][5] != 1
                            boolFinished = false;
                        end
                    end
                    println("-------------",boolFinished," ",keepIter," ",[treeList[l][5] for l in 1:length(treeList)],[treeList[l][3] for l in 1:length(treeList)],"-------------");
                    if (boolFinished) || (!(keepIter))
                        println("break");
                        break
                    end
                    openNodes = [(treeList[l][3],l) for l in 1:length(treeList) if treeList[l][5] == -1];
                    if openNodes != []
                        selectNode = sort(openNodes, by = x -> x[1])[1][2];
                        println("On core: ",p," processing node: ",selectNode);
                        treeList[selectNode][5] = 0;
                        mpSolveInfo = remotecall_fetch(solveMP_para_Share,p,[selectNode,tcore,xcore,weigthCore,ubCost,
                            tbest,xbest,noTh,wpList]);
                        # update the cutList with the added cuts and two new nodes
                        # update the cutSet
                        if mpSolveInfo[6] < ubCost
                            ubCost = mpSolveInfo[6];
                            tbest = mpSolveInfo[4];
                            xbest = mpSolveInfo[5];
                        end
                        append!(tcoreList,mpSolveInfo[8]);
                        append!(xcoreList,mpSolveInfo[9]);
                        append!(ycoreList,mpSolveInfo[10]);
                        append!(ubcoreList,mpSolveInfo[11]);
                        # compare the current lb with the current best lb
                        lbOverAll = minimum([treeList[l][3] for l in 1:length(treeList) if treeList[l][5] != 1]);
                        mpStatus = mpSolveInfo[1];
                        mpObj = mpSolveInfo[2];
                        if mpStatus == :Optimal
                            if mpObj < lbOverAll
                                lbOverAll = mpObj;
                            end
                        end
                        if (ubCost - lbOverAll)/ubCost < ϵ
                            keepIter = false;
                        else
                            if (mpSolveInfo[2] < ubCost)&(mpSolveInfo[1] == :Optimal)
                                # branch the current node
                                GList = mpSolveInfo[3];
                                cutSet = mpSolveInfo[7];
                                GCurrent = GList[length(GList)];
                                GFrac = Dict();
                                lbCost = mpSolveInfo[2];
                                for i in pData.II
                                    GFraciList = [disData[ω].H for ω in Ω if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)];
                                    if GFraciList != []
                                        Gbegin = findfirst(HShare,minimum(GFraciList));
                                        Gend = findfirst(HShare,maximum(GFraciList));
                                        GFrac[i] = [Gbegin,Gend];
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
                                # select the largest GFrac coverage
                                largestGFrac = -Inf;
                                lGFracInd = -1;
                                for i in pData.II
                                    if GFrac[i] != []
                                        if GFrac[i][2] - GFrac[i][1] > largestGFrac
                                            largestGFrac = GFrac[i][2] - GFrac[i][1];
                                            lGFracInd = i;
                                        end
                                    end
                                end
                                locBreak = Int64(floor((GFrac[lGFracInd][1]*2/3 + GFrac[lGFracInd][2]*1/3)));
                                divSet1,divDet1,divSet2,divDet2 = breakDiv(pData,disData,H,treeList[selectNode][1],treeList[selectNode][2],lGFracInd,locBreak,distanceDict);
                                lGFracInd1 = -1;
                                largest1 = -Inf;
                                fracTopLG = -1;
                                fracBotLG = -1;
                                for i in pData.II
                                    if (!(i in allSucc[lGFracInd]))&(!(lGFracInd in allSucc[i]))&(i != lGFracInd)
                                        if GFrac[i] != []
                                            fracTop = min(GFrac[i][2],maximum([divSet1[i][par].endH for par in 1:length(divSet1[i]) if divDet1[i][par] == 0]));
                                            fracBot = max(GFrac[i][1],minimum([divSet1[i][par].startH for par in 1:length(divSet1[i]) if divDet1[i][par] == 0]));
                                            #println(i," ",fracTop," ",fracBot);
                                            if fracTop - fracBot > largest1
                                                largest1 = fracTop - fracBot;
                                                lGFracInd1 = i;
                                                fracTopLG = fracTop;
                                                fracBotLG = fracBot;
                                            end
                                        end
                                    end
                                end
                                newPartition1 = [(lGFracInd1,fracBotLG,fracTopLG)];
                                divSet1,divDet1 = splitPar3(divSet1,divDet1,newPartition1);
                                println(newPartition1);
                                push!(treeList,[divSet1,divDet1,lbCost,deepcopy(cutSet),-1]);

                                lGFracInd2 = -1;
                                largest2 = -Inf;
                                fracTopLG = -1;
                                fracBotLG = -1;
                                for i in pData.II
                                    if (!(i in allSucc[lGFracInd]))&(!(lGFracInd in allSucc[i]))&(i != lGFracInd)
                                        if GFrac[i] != []
                                            fracTop = min(GFrac[i][2],maximum([divSet2[i][par].endH for par in 1:length(divSet2[i]) if divDet2[i][par] == 0]));
                                            fracBot = max(GFrac[i][1],minimum([divSet2[i][par].startH for par in 1:length(divSet2[i]) if divDet2[i][par] == 0]));
                                            if fracTop - fracBot > largest2
                                                largest2 = fracTop - fracBot;
                                                lGFracInd2 = i;
                                                fracTopLG = fracTop;
                                                fracBotLG = fracBot;
                                            end
                                        end
                                    end
                                end
                                newPartition2 = [(lGFracInd2,fracBotLG,fracTopLG)];
                                divSet2,divDet2 = splitPar3(divSet2,divDet2,newPartition2);
                                println(newPartition2);
                                push!(treeList,[divSet2,divDet2,lbCost,deepcopy(cutSet),-1]);

                                treeList[selectNode][5] = 1;
                            else
                                treeList[selectNode][5] = 1;
                            end
                        end
                    else
                        # sleep if there is no current available nodes
                        println(p);
                        remotecall_fetch(sleep, p, 30);
                    end
                end
            end
        end
    end
    return tbest,xbest,ubCost,lbOverAll;
end

function partSolve_BB_para(pData,disData,Ω,sN,MM,noThreads,ϵ = 1e-2)
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
    lDictShare = SharedArray{Float64,1}((length(pData.II)));
    for i in 1:length(pData.II)
        lDict[pData.II[i]] += disData[length(Ω)].H;
        lDictShare[i] = lDict[pData.II[i]];
    end
    Tmax1 = lDict[0];

    allSucc = findSuccAll(pData);
    # the first two indices are float, needs to be convert to integer when being used
    distanceShare = SharedArray{Float64,2}((sum(length(allSucc[iKey]) for iKey in keys(allSucc)),3));
    iNo = 0;
    for i in pData.II
        for j in allSucc[i]
            iNo += 1;
            distanceShare[iNo,1] = i;
            distanceShare[iNo,2] = j;
            distanceShare[iNo,3] = detCal(pData,i,j);
        end
    end

    # make the data SharedArray
    p0BShare = SharedArray{Float64,1}(2);
    p0BShare[1] = pData.p0;
    p0BShare[2] = pData.B;
    IIShare = SharedArray{Int,1}(length(pData.II));
    DShare = SharedArray{Int,1}(length(pData.II));
    for i in 1:length(pData.II)
        IIShare[i] = pData.II[i];
        DShare[i] = pData.D[pData.II[i]];
    end
    KShare = SharedArray{Int,2}((length(pData.K),2));
    for k in 1:length(pData.K)
        KShare[k,1] = pData.K[k][1];
        KShare[k,2] = pData.K[k][2];
    end
    bShare = SharedArray{Float64,2}((length(pData.II),maximum([length(values(pData.Ji[i])) for i in pData.II])));
    effShare = SharedArray{Float64,2}((length(pData.II),maximum([length(values(pData.Ji[i])) for i in pData.II])));
    for i in 1:length(pData.II)
        for j in 1:length(pData.Ji[pData.II[i]])
            bShare[i,j] = pData.b[pData.II[i]][pData.Ji[pData.II[i]][j]];
            effShare[i,j] = pData.eff[pData.II[i]][pData.Ji[pData.II[i]][j]];
        end
    end

    HOriShare = SharedArray{Float64,1}(length(Ω));
    for ω in 1:length(Ω)
        HOriShare[ω] = disData[Ω[ω]].H;
    end
    disdShare = SharedArray{Float64,2}((length(pData.II),length(Ω)));
    disPrShare = SharedArray{Float64,1}((length(Ω)));
    for ω in 1:length(Ω)
        for i in 1:length(pData.II)
            if pData.II[i] != 0
                disdShare[i,ω] = disData[Ω[ω]].d[pData.II[i]];
            else
                disdShare[i,ω] = 0;
            end
        end
        disPrShare[ω] = disData[Ω[ω]].prDis;
    end
    # recover the information and make them everywhere
    @everywhere global pData = recoverpInfo(IIShare,DShare,effShare,bShare,KShare,p0BShare[1],p0BShare[2]);
    @everywhere global disData; global H; global Ω = recoverdInfo(IIShare,HOriShare,disdShare,disPrShare);
    @everywhere global distanceDict; global allSucc = recoverdis(distanceShare);
    @everywhere global lDict = recoverlDict(IIShare,lDictShare);

    # start with an upper bound based on the smaller stochastic solution
    ubextList,tHList,ubInc,tbest,xbest,θbest,textList,xextList = iniPart(pData,disData,Ω,sN,MM,1,noThreads);
    lbCost = -Inf;
    lbCostList = [];
    global ubCost = ubInc;

    IJPair = [(i,j) for i in pData.II for j in pData.Ji[i]];
    textShare = SharedArray{Float64,2}((length(pData.II),length(textList)));
    for i in pData.II
        for it in 1:length(textList)
            textShare[i,it] = textList[it][i];
        end
    end
    xextShare = SharedArray{Float64,2}((length(IJPair),length(xextList)));
    counter = 1;
    for i in pData.II
        for j in pData.Ji
            for itj in 1:length(xextList)
                xextShare[counter,itj] = xextList[itj][IJPair[counter]];
            end
            counter += 1;
        end
    end
    ubextShare = SharedArray{Float64,1}(length(ubextList));
    for itu in 1:length(ubextList)
        ubextShare[itu] = ubextList[itu];
    end
    @everywhere global tcoreList; global xcoreList; global ubextList = recoverCoreList(IIShare,textShare,xextShare,ubextShare);

    brInfo = precludeRelNew(pData,H,ubCost);

    # initialize cutSet and divSet
    HΩ = 1:(length(H) - 2);
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

    IPPair = [(i,par) for i in pData.II for par in 1:length(divSet[i])];
    πSet = SharedArray{Float64,3}((0,length(pData.II),length(Ω)));
    λSet = SharedArray{Float64,3}((0,length(IJPair),length(Ω)));
    γSet = SharedArray{Float64,3}((0,length(IPPair),length(Ω)));
    vSet = SharedArray{Float64,2}((0,length(Ω)));
    thatSet = SharedArray{Float64,2}((0,length(pData.II)));
    xhatSet = SharedArray{Float64,2}((0,length(IJPair)));
    yhatSet = SharedArray{Int,2}((0,length(IPPair)));
    divInfoShare = SharedArray{Int,2}((length(IPPair),5));
    for ip in 1:length(IPPair)
        divInfoShare[ip,1] = IPPair[ip][1];
        divInfoShare[ip,2] = IPPair[ip][2];
        divInfoShare[ip,3] = divSet[divInfoShare[ip,1]][divInfoShare[ip,2]].startH;
        divInfoShare[ip,4] = divSet[divInfoShare[ip,1]][divInfoShare[ip,2]].endH;
        divInfoShare[ip,5] = divDet[divInfoShare[ip,1]][divInfoShare[ip,2]];
    end

    # set up a tree list
    @everywhere global treeList = [];
    @everywhere global cutList = [];
    @everywhere push!(treeList,[lbCost,[],-1]); # the empty set is the list of predecessors of the current node
    @everywhere push!(cutList,[πSet,λSet,γSet,vSet,thatSet,xhatSet,yhatSet,divInfoShare]);

    global batchNo = 5;
    global lbOverAll = -Inf;
    # transfer the data back to everywhere
    tbest,xbest,ubCost,lbOverAll = runPara_Share(treeList,cutList,tcoreList,xcoreList,ubCost,tbest,xbest,batchNo);

    # need a cut selection process within the callback
    return tbest,xbest,ubCost,lbOverAll;
end
