function recoverpInfo(II,D,eff,b,K,p0,B)
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

function recoverdInfo(II,HOriShare,disdShare,disPrShare,Tmax)
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
    HRev = Dict();
    for hIter in keys(H)
        HRev[H[hIter]] = hIter;
    end
    return disData,H,Ω,HRev;
end

function recoverdis(IIShare,distanceShare)
    allSucc = Dict();
    distanceDict = Dict();
    for i in IIShare
        allSucc[i] = [];
    end
    for i in 1:size(distanceShare)[1]
        push!(allSucc[Int(distanceShare[i,1])],Int(distanceShare[i,2]));
        distanceDict[Int(distanceShare[i,1]),Int(distanceShare[i,2])] = distanceShare[i,3];
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

function convertDiv(divSet,divDet)
    IPPair = [(i,par) for i in pData.II for par in 1:length(divSet[i])];
    divInfoShare = SharedArray{Int,2}((length(IPPair),5));
    for ip in 1:length(IPPair)
        divInfoShare[ip,1] = IPPair[ip][1];
        divInfoShare[ip,2] = IPPair[ip][2];
        divInfoShare[ip,3] = divSet[divInfoShare[ip,1]][divInfoShare[ip,2]].startH;
        divInfoShare[ip,4] = divSet[divInfoShare[ip,1]][divInfoShare[ip,2]].endH;
        divInfoShare[ip,5] = divDet[divInfoShare[ip,1]][divInfoShare[ip,2]];
    end
    return divInfoShare;
end

function recoverDiv(divInfoShare)
    divSet = Dict();
    divDet = Dict();
    for di in 1:size(divInfoShare)[1]
        if !(divInfoShare[di,1] in keys(divSet))
            divSet[divInfoShare[di,1]] = [partType(divInfoShare[di,3],divInfoShare[di,4])];
            divDet[divInfoShare[di,1]] = [divInfoShare[di,5]];
        else
            push!(divSet[divInfoShare[di,1]],partType(divInfoShare[di,3],divInfoShare[di,4]));
            push!(divDet[divInfoShare[di,1]],divInfoShare[di,5]);
        end
    end
    return divSet,divDet;
end

function recoverCoreList(II,IJPair,textShare,xextShare,ubextShare)
    textList = [];
    xextList = [];
    ubextList = [];
    for it in 1:size(textShare)[2]
        tDict = Dict();
        for i in 1:length(II)
            tDict[II[i]] = textShare[i,it];
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

function runRecover(IIShare,DShare,effShare,bShare,KShare,p0BShare,HOriShare,disdShare,disPrShare,distanceShare,lDictShare)
    # recover the information and make them everywhere
    global pData = recoverpInfo(IIShare,DShare,effShare,bShare,KShare,p0BShare[1],p0BShare[2]);
    global (disData, H, Ω, HRev) = recoverdInfo(IIShare,HOriShare,disdShare,disPrShare,p0BShare[3]);
    global (allSucc,distanceDict) = recoverdis(IIShare,distanceShare);
    global lDict = recoverlDict(IIShare,lDictShare);
end

function subPara1(pData,disData,Ω,tbest,xbest,ybest,divSet,H,lDict,wp = CachingPool(workers()))
    θList = pmap(wp,ω -> sub_divT(pData,disData[ω],ω,tbest,xbest,ybest,divSet,H,lDict),Ω);
    return θList;
end

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
            ilb = divSet[i][minimum(ibSet)].startH;
            if (tSol[i] < H[ilb])|(tSol[i] >= H[iub])
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
    # input: [cutData,cutCurrent,tcoreList,xcoreList,ubcoreList,ubCost,tbest,xbest,noTh,wpList]
    divSet,divDet = recoverDiv(data[3]);
    divData = data[2];
    cutSet = data[1];           # historical cuts
    IJPair = [(i,j) for i in pData.II for j in pData.Ji[i]];
    IPPair = [(i,par) for i in pData.II for par in 1:length(divSet[i])];
    tcoreList,xcoreList,ubcoreList = recoverCoreList(pData.II,IJPair,data[4],data[5],data[6]);
    ubCost = data[7];
    tbest = data[8];
    xbest = data[9];
    noTh = data[10];
    wp = CachingPool(data[11]);

    Tmax1 =lDict[0];
    GList = [];
    tcoreNew = [];
    xcoreNew = [];
    ubcoreNew = [];
    cutSetNew = [];

    function partBenders(cb)
        currentLB = MathProgBase.cbgetbestbound(cb);
        println("lazy,$(currentLB)");
        if currentLB <= minimum(ubCostList)
            # the callback function
            that = SharedArray{Float64,1}((length(pData.II)));
            tdict = Dict();
            xhat = SharedArray{Float64,1}((length(IJPair)));
            xdict = Dict();
            θhat = SharedArray{Float64,1}((length(Ω)));
            yhat = SharedArray{Float64,1}((length(IPPair)));
            ydict = Dict();
            # obtain the solution at the current node
            for i in pData.II
                that[findfirst(pData.II,i)] = getvalue(t[i]);
                tdict[i] = getvalue(t[i]);
                for j in pData.Ji[i]
                    xhat[findfirst(IJPair,(i,j))] = getvalue(x[i,j]);
                    xdict[i,j] = getvalue(x[i,j]);
                end
                for par in 1:length(divSet[i])
                    yhat[findfirst(IPPair,(i,par))] = round(getvalue(y[i,par]));
                    ydict[i,par] = round(getvalue(y[i,par]));
                end
            end
            for ω in 1:length(Ω)
                θhat[ω] = getvalue(θ[Ω[ω]]);
            end
            push!(tcoreList,tdict);
            push!(xcoreList,xdict);
            push!(tcoreNew,that);
            push!(xcoreNew,xhat);

            # generate cuts
            θInt = Dict();
            ubCost = minimum(ubCostList);
            ubTemp,θInt = ubCalP(pData,disData,Ω,xdict,tdict,Tmax1,1,wp);
            push!(ubcoreList,ubTemp);
            push!(ubcoreNew,ubTemp);

            if ubCost > ubTemp
                for i in pData.II
                    tbest[i] = that[findfirst(pData.II,i)];
                    for j in pData.Ji[i]
                        xbest[i,j] = xhat[findfirst(IJPair,(i,j))];
                    end
                end
            end
            push!(ubCostList,ubTemp);

            #dataList = pmap(ω -> sub_divT(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict), Ω);
            # obtain the cores
            tcore,xcore,ycore = avgCore(pData,divSet,tcoreList,xcoreList,ycoreList);
            # here is the issue, pack it in a function prevent separating it
            dataList = subPara(pData,disData,Ω,tdict,xdict,ydict,divSet,H,lDict,tcore,xcore,ycore,wp);
            cutScen = [ω for ω in Ω if dataList[ω][4] - θhat[findfirst(Ω,ω)] > 1e-4*θhat[findfirst(Ω,ω)]];
            πSet = SharedArray{Float64,2}((length(pData.II),length(cutScen)));
            λSet = SharedArray{Float64,2}((length(IJPair),length(cutScen)));
            γSet = SharedArray{Float64,2}((length(IPPair),length(cutScen)));
            vSet = SharedArray{Float64,1}((length(cutScen)));
            for ωi in 1:length(cutScen)
                ω = cutScen[ωi];
                vSet[ωi] = dataList[ω][4];
                for i in pData.II
                    πSet[findfirst(pData.II,i),ωi] = dataList[ω][1][i];
                    for j in pData.Ji[i]
                        λSet[findfirst(IJPair,(i,j)),ωi] = dataList[ω][2][i,j];
                    end
                    for par in 1:length(divSet[i])
                        γSet[findfirst(IPPair,(i,par)),ωi] = dataList[ω][3][i,par];
                    end
                end
                @lazyconstraint(cb, θ[ω] >= vSet[ωi] + sum(πSet[findfirst(pData.II,i),ωi]*(t[i] - tdict[i]) for i in pData.II) +
                    sum(sum(λSet[findfirst(IJPair,(i,j)),ωi]*(x[i,j] - xdict[i,j]) for j in pData.Ji[i]) for i in pData.II) +
                    sum(sum(γSet[findfirst(IPPair,(i,par)),ωi]*(y[i,par] - ydict[i,par]) for par in 1:length(divSet[i])) for i in pData.II));
            end
            newCuts = [cutScen,πSet,λSet,γSet,vSet,that,xhat,yhat];
            #push!(cutSet,[[that,xhat,yhat,divSet],cutDual]);
            push!(cutSetNew,newCuts);
            GCurrent = [dataList[ω][5] for ω in Ω];
            push!(GList,GCurrent);
        else
            return JuMP.StopTheSolver;
        end
    end

    # correct all ycoreList
    ubCostList = [ubCost];
    ycoreList = Dict();
    for ll in 1:length(tcoreList)
        ycoreList[ll] = Dict();
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
        xfeas = xcoreList[tcoreInd];
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
    # nc is the node number
    # npoint is the solution index
    # cutSet[nc] = [npoint,[],πSet,λSet,γSet,vSet,thatSet,xhatSet,yhatSet,divInfoShare]
    for nc in 1:length(cutSet)
        divSetPrev,divDetPrev = recoverDiv(divData[nc]);
        IPPairPrev = [(i,par) for i in pData.II for par in 1:length(divSetPrev[i])];
        for npoint in 1:length(cutSet[nc])
            for ωi in cutSet[nc][npoint][1]
                ω = cutSet[nc][npoint][1][ωi];
                vk = cutSet[nc][npoint][5][ωi];
                πk = cutSet[nc][npoint][2][:,ωi];
                λk = cutSet[nc][npoint][3][:,ωi];
                γk = cutSet[nc][npoint][4][:,ωi];
                @constraint(mp, θ[ω] >= vk + sum(πk[findfirst(pData.II,i)]*(mp[:t][i] - cutSet[nc][npoint][6][findfirst(pData.II,i)]) +
                    sum(λk[findfirst(IJPair,(i,j))]*(mp[:x][i,j] - cutSet[nc][npoint][7][findfirst(IJPair,(i,j))]) for j in pData.Ji[i]) +
                    sum(γk[findfirst(IPPairPrev,(i,par))]*(sum(mp[:y][i,parNew] for parNew in 1:length(divSet[i]) if revPar(divSetPrev[i],divSet[i][parNew]) == par) - cutSet[nc][npoint][8][findfirst(IPPairPrev,(i,par))])
                    for par in 1:length(divSetPrev[i])) for i in pData.II));
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

    if mpStatus == :Optimal
        returnNo = mpObj;
        # branch
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
        divSet1,divDet1,divSet2,divDet2 = breakDiv(pData,disData,H,divSet,divDet,lGFracInd,locBreak,distanceDict);
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
        divShare1 = convertDiv(divSet1,divDet1);

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
        divShare2 = convertDiv(divSet2,divDet2);

        returnSet = [divShare1,divShare2];
    else
        returnNo = -Inf;
        returnSet = [];
    end

# mpStatus,mpObj,GList,tbest,xbest,minimum(ubCostList),cutSet,tcoreNew,xcoreNew,ycoreNew,ubcoreNew;
    return returnNo,cutSetNew,returnSet,tbest,xbest,minimum(ubCostList);
end

function runPara_Share(treeList,cutList,tcoreShare,xcoreShare,ubcoreShare,ubCost,tbest,xbest,batchNo)
    npList = workers()[1:batchNo];
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
                        if treeList[l][3] != 1
                            boolFinished = false;
                        end
                    end
                    println("-------------",boolFinished," ",keepIter," ",[treeList[l][3] for l in 1:length(treeList)],[treeList[l][1] for l in 1:length(treeList)],"-------------");
                    if (boolFinished) || (!(keepIter))
                        println("break");
                        break
                    end
                    openNodes = [(treeList[l][1],l) for l in 1:length(treeList) if treeList[l][3] == -1];
                    if openNodes != []
                        selectNode = sort(openNodes, by = x -> x[1])[1][2];
                        println("On core: ",p," processing node: ",selectNode);
                        treeList[selectNode][3] = 0;
                        cutData = cutList[treeList[selectNode][2]];
                        divData = [treeList[id][4] for id in treeList[selectNode][2]];
                        mpSolveInfo = remotecall_fetch(solveMP_para_Share,p,[cutData,divData,treeList[selectNode][4],tcoreShare,xcoreShare,ubcoreShare,ubCost,tbest,xbest,noTh,wpList]);
                        # update the cutList with the added cuts and two new nodes
                        # update the cutSet
                        # return returnNo,cutSet,returnSet,tbest,xbest,minimum(ubCostList)
                        lbOverAll = minimum([treeList[l][1] for l in 1:length(treeList) if treeList[l][3] != 1]);
                        if mpSolveInfo[1] > -Inf
                            if mpSolveInfo[6] < ubCost
                                ubCost = mpSolveInfo[6];
                                tbest = mpSolveInfo[4];
                                xbest = mpSolveInfo[5];
                            end
                            if mpSolveInfo[1] < lbOverAll
                                lbOverAll = mpSolveInfo[1];
                            end
                            if (ubCost - lbOverAll)/ubCost < ϵ
                                keepIter = false;
                            else
                                # update the current node cut
                                cutList[selectNode] = mpSolveInfo[2];
                                # push the branched nodes
                                if (mpSolveInfo[1] < ubCost)
                                    ancestorTemp = deepcopy(treeList[selectNode][2]);
                                    push!(ancestorTemp,selectNode);
                                    for newN in 1:length(mpSolveInfo[3])
                                        push!(treeList,[mpSolveInfo[1],ancestorTemp,-1,mpSolveInfo[3][newN]]);
                                        push!(cutList,[]);
                                    end
                                end
                            end
                        end
                        treeList[selectNode][3] = 1;
                    else
                        println("**************Worker $(p) waiting for open nodes.**************");
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
    p0BShare = SharedArray{Float64,1}(3);
    p0BShare[1] = pData.p0;
    p0BShare[2] = pData.B;
    p0BShare[3] = Tmax;
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

    for i in procs()
       remotecall_fetch(runRecover,i,IIShare,DShare,effShare,bShare,KShare,p0BShare,HOriShare,disdShare,disPrShare,distanceShare,lDictShare);
    end

    # start with an upper bound based on the smaller stochastic solution
    ubextList,tHList,ubInc,tbest,xbest,θbest,textList,xextList = iniPart(pData,disData,Ω,sN,MM,1,noThreads);
    lbCost = -Inf;
    lbCostList = [];
    global ubCost = ubInc;

    IJPair = [(i,j) for i in pData.II for j in pData.Ji[i]];
    textShare = SharedArray{Float64,2}((length(pData.II),length(textList)));
    for i in 1:length(pData.II)
        for it in 1:length(textList)
            textShare[i,it] = textList[it][pData.II[i]];
        end
    end
    xextShare = SharedArray{Float64,2}((length(IJPair),length(xextList)));
    counter = 1;
    for i in pData.II
        for j in pData.Ji[i]
            for itj in 1:length(xextList)
                xextShare[counter,itj] = xextList[itj][i,j];
            end
            counter += 1;
        end
    end
    ubextShare = SharedArray{Float64,1}(length(ubextList));
    for itu in 1:length(ubextList)
        ubextShare[itu] = ubextList[itu];
    end

    brInfo = precludeRelNew(pData,H,ubCost);

    # initialize cutSet and divSet
    HΩ = 1:length(H) - 2;
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

    divInfoShare = convertDiv(divSet,divDet);

    # set up a tree list
    global treeList = [];
    global cutList = [];
    push!(treeList,[lbCost,[],-1,divInfoShare]); # the empty set is the list of predecessors of the current node
    npoint = 0;
    push!(cutList,[]);

    global batchNo = 5;
    global lbOverAll = -Inf;
    # transfer the data back to everywhere
    tbest,xbest,ubCost,lbOverAll = runPara_Share(treeList,cutList,textShare,xextShare,ubextShare,ubCost,tbest,xbest,batchNo);

    # need a cut selection process within the callback
    return tbest,xbest,ubCost,lbOverAll;
end
