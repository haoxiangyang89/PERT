# B&B process with both MW and UB

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
    #divSet,divDet = recoverDiv(data[3]);
    divSet,divDet = data[3];
    divData = data[2];
    cutSet = data[1];           # historical cuts
    IJPair = [(i,j) for i in pData.II for j in pData.Ji[i]];
    IPPair = [(i,par) for i in pData.II for par in 1:length(divSet[i])];
    #tcoreList,xcoreList,ubcoreList = recoverCoreList(pData.II,IJPair,data[4],data[5],data[6]);
    tcoreList = data[4];
    xcoreList = data[5];
    ubcoreList = data[6];
    ubCost = data[7];
    tbest = data[8];
    xbest = data[9];
    noTh = data[10];
    wp = CachingPool(data[11]);
    nSplit = data[12];
    roundLimit = data[13];
    cutSelOpt = data[14];
    Ω = 1:length(disData);

    Tmax1 =lDict[0];
    GList = [];
    tcoreNew = [];
    xcoreNew = [];
    ycoreNew = [];
    ubcoreNew = [];
    cutSetNew = [];
    tError = [];
    xError = [];
    yError = [];
    ErrorωList = [];
    tcoreError = [];
    xcoreError = [];
    ycoreError = [];
    tUnbounded = [];
    xUnbounded = [];
    yUnbounded = [];
    UnboundedωList = [];
    tcoreUnbounded = [];
    xcoreUnbounded = [];
    ycoreUnbounded = [];

    function partBenders(cb)
        currentLB = MathProgBase.cbgetbestbound(cb);
        println("lazy,$(currentLB)");
        if currentLB <= minimum(ubCostList)
            # the callback function
            that = Dict();
            xhat = Dict();
            yhat = Dict();
            θhat = Dict();
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
            for ω in 1:length(Ω)
                θhat[ω] = getvalue(θ[Ω[ω]]);
            end
            push!(tcoreList,that);
            push!(xcoreList,xhat);
            push!(ycoreList,yhat);
            push!(tcoreNew,that);
            push!(xcoreNew,xhat);
            push!(ycoreNew,yhat);

            # generate cuts
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

            # obtain the cores
            tcore,xcore,ycore = avgCore(pData,divSet,tcoreList,xcoreList,ycoreList);
            # here is the issue, pack it in a function prevent separating it
            dataList = subPara(pData,disData,Ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore,wp);
            cutScen = [];
            errorInd = false;
            for ω in Ω
                if length(dataList[ω]) == 3
                    push!(tError,that);
                    push!(xError,xhat);
                    push!(yError,yhat);
                    push!(ErrorωList,ω);
                    push!(tcoreError,tcore);
                    push!(xcoreError,xcore);
                    push!(ycoreError,ycore);
                    errorInd = true;
                else
                    if (dataList[ω][4] - θhat[findfirst(Ω,ω)] > 1e-4*θhat[findfirst(Ω,ω)])
                        push!(cutScen,ω);
                    end
                    if length(dataList[ω]) == 6
                        push!(tUnbounded,that);
                        push!(xUnbounded,xhat);
                        push!(yUnbounded,yhat);
                        push!(UnboundedωList,ω);
                        push!(tcoreUnbounded,tcore);
                        push!(xcoreUnbounded,xcore);
                        push!(ycoreUnbounded,ycore);
                    end
                end
            end
            if !(errorInd)
                GCurrent = [dataList[ω][5] for ω in Ω];
                push!(GList,GCurrent);
            end
            πSet = zeros(length(pData.II),length(cutScen));
            λSet = zeros(length(IJPair),length(cutScen));
            γSet = zeros(length(IPPair),length(cutScen));
            πSet1 = zeros(length(pData.II),length(cutScen));
            λSet1 = zeros(length(IJPair),length(cutScen));
            γSet1 = zeros(length(IPPair),length(cutScen));
            vSet = zeros(length(cutScen));
            for ωi in 1:length(cutScen)
                ω = cutScen[ωi];
                vSet[ωi] = dataList[ω][4];
                for i in pData.II
                    vSet[ωi] -= dataList[ω][1][i]*that[i];
                    if abs(dataList[ω][1][i]) >= 1e-7
                        πSet[findfirst(pData.II,i),ωi] = dataList[ω][1][i];
                    else
                        πSet[findfirst(pData.II,i),ωi] = 0;
                        if dataList[ω][1][i] < 0
                            vSet[ωi] += dataList[ω][1][i];
                        end
                    end
                    for j in pData.Ji[i]
                        vSet[ωi] -= dataList[ω][2][i,j]*xhat[i,j];
                        if abs(dataList[ω][2][i,j]) >= 1e-5
                            λSet[findfirst(IJPair,(i,j)),ωi] = dataList[ω][2][i,j];
                        else
                            λSet[findfirst(IJPair,(i,j)),ωi] = 0;
                            if dataList[ω][2][i,j] < 0
                                vSet[ωi] += dataList[ω][2][i,j];
                            end
                        end
                    end
                    for par in 1:length(divSet[i])
                        vSet[ωi] -= dataList[ω][3][i,par]*yhat[i,par];
                        if abs(dataList[ω][3][i,par]) >= 1e-5
                            γSet[findfirst(IPPair,(i,par)),ωi] = dataList[ω][3][i,par];
                        else
                            γSet[findfirst(IPPair,(i,par)),ωi] = 0;
                            if dataList[ω][3][i,par] < 0
                                vSet[ωi] += dataList[ω][3][i,par];
                            end
                        end
                    end
                end
                @lazyconstraint(cb, θ[ω] >= vSet[ωi] + sum(πSet[findfirst(pData.II,i),ωi]*t[i] for i in pData.II) +
                    sum(sum(λSet[findfirst(IJPair,(i,j)),ωi]*x[i,j] for j in pData.Ji[i]) for i in pData.II) +
                    sum(sum(γSet[findfirst(IPPair,(i,par)),ωi]*y[i,par] for par in 1:length(divSet[i])) for i in pData.II));
            end
            newCuts = [cutScen,πSet,λSet,γSet,vSet];
            #push!(cutSet,[[that,xhat,yhat,divSet],cutDual]);
            push!(cutSetNew,newCuts);
        else
            return JuMP.StopTheSolver;
        end
    end

    # correct all ycoreList
    ubCostList = [ubCost];
    ycoreList = [];
    for ll in 1:length(tcoreList)
        yTemp = Dict();
        for i in pData.II
            for par in 1:length(divSet[i])
                if (tcoreList[ll][i] >= H[divSet[i][par].startH])&(tcoreList[ll][i] < H[divSet[i][par].endH])
                    yTemp[i,par] = 1;
                else
                    yTemp[i,par] = 0;
                end
            end
        end
        push!(ycoreList,yTemp);
    end

    # move the createMaster_Callback here
    mp = Model(solver = GurobiSolver(GUROBI_ENV,IntFeasTol = 1e-8, FeasibilityTol = 1e-8, Threads = noTh, Cutoff = ubCost, TimeLimit = roundLimit));
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
    # else
    #     # solve the simple problem to identify a feasible solution
    #     # combine all 0 partitions to be a single partition, solve the master problem on this simple partition
    #     tfeas,xfeas = combinePart(pData,disData,Ω,divSet,divDet,H,tcoreList,xcoreList,ycoreList,wp,noTh,ubCost,Tmax1);
    #     yfeas = Dict();
    #     for i in pData.II
    #         setvalue(t[i], tfeas[i]);
    #         for j in pData.Ji[i]
    #             setvalue(x[i,j],xfeas[i,j]);
    #         end
    #         for par in 1:length(divSet[i])
    #             if (tfeas[i] >= H[divSet[i][par].startH])&(tfeas[i] < H[divSet[i][par].endH])
    #                 yfeas[i,par] = 1;
    #             elseif (abs(tfeas[i] - H[length(H) - 1]) < 1e-4)&(divSet[i][par].endH == length(H) - 1)
    #                 yfeas[i,par] = 1;
    #             else
    #                 yfeas[i,par] = 0;
    #             end
    #             setvalue(y[i,par],yfeas[i,par]);
    #         end
    #     end
    #     θfeas = subPara1(pData,disData,Ω,tfeas,xfeas,yfeas,divSet,H,lDict,wp);
    #     for ω in Ω
    #         setvalue(θ[ω],θfeas[ω]);
    #     end
    end

    tCurrent = Dict();
    xCurrent = Dict();
    θCurrent = Dict();
    yCurrent = Dict();
    GCurrent = Dict();

    # add the cut
    # cutInfo = 2 dimensional vector, first dimention record the primal solution,
    # second dimension record the dual solution for every scenario
    # nc is the node number
    # npoint is the solution index
    # cutSet[nc] = [npoint,[],πSet,λSet,γSet,vSet,thatSet,xhatSet,yhatSet,divInfoShare]
    for nc in 1:length(cutSet)
        #divSetPrev,divDetPrev = recoverDiv(divData[nc]);
        divSetPrev,divDetPrev = divData[nc];
        revDict = Dict();
        for i in pData.II
            revDict[i] = Dict();
            for par in 1:length(divSetPrev[i])
                revDict[i][par] = [];
                for parNew in 1:length(divSet[i])
                    if (divSet[i][parNew].startH >= divSetPrev[i][par].startH)&(divSet[i][parNew].endH <= divSetPrev[i][par].endH)
                        push!(revDict[i][par],parNew);
                    end
                end
            end
        end
        IPPairPrev = [(i,par) for i in pData.II for par in 1:length(divSetPrev[i])];
        for npoint in 1:length(cutSet[nc])
            for ωi in 1:length(cutSet[nc][npoint][1])
                ω = cutSet[nc][npoint][1][ωi];
                vk = cutSet[nc][npoint][5][ωi];
                πk = cutSet[nc][npoint][2][:,ωi];
                λk = cutSet[nc][npoint][3][:,ωi];
                γk = cutSet[nc][npoint][4][:,ωi];
                @constraint(mp, θ[ω] >= vk + sum(πk[findfirst(pData.II,i)]*mp[:t][i] +
                    sum(λk[findfirst(IJPair,(i,j))]*mp[:x][i,j] for j in pData.Ji[i]) +
                    sum(γk[findfirst(IPPairPrev,(i,par))]*(sum(mp[:y][i,parNew] for parNew in revDict[i][par]))
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

    if (mpStatus == :Optimal)
        returnNo = mpObj;
        for i in pData.II
            tCurrent[i] = getvalue(mp[:t][i]);
            for j in pData.Ji[i]
                xCurrent[i,j] = getvalue(mp[:x][i,j]);
            end
            for par in 1:length(divSet[i])
                yCurrent[i,par] = getvalue(mp[:y][i,par]);
            end
        end
        subInfo = pmap(wp,ω -> sub_divT(pData,disData[ω],ω,tCurrent,xCurrent,yCurrent,divSet,H,lDict,2),Ω);
        GCurrent = [subInfo[ω][2] for ω in Ω];
        θCurrent = [subInfo[ω][1] for ω in Ω];
        ubCurrent,θIntCurrent = ubCalP(pData,disData,Ω,xCurrent,tCurrent,Tmax1,1,wp);
        # branch
        θDiff = [θIntCurrent[ω] - θCurrent[ω] for ω in Ω];
        θDiffPerm = sortperm(θDiff,rev = true);
        locBreak = θDiffPerm[1];
        locBreakH = HRev[disData[locBreak].H];
        θgain = 0;
        lGFracInd = -1;
        for i in pData.II
            θgainTemp = 0;
            if i != 0
                if tCurrent[i] > H[locBreakH] + 1e-6
                    if (GCurrent[locBreak][i] < 1 - 1e-6)&(GCurrent[locBreak][i] > 1e-6)
                        θgainTemp += (1 - GCurrent[locBreak][i])*disData[locBreak].d[i];
                    end
                elseif tCurrent[i] < H[locBreakH] - 1e-6
                    if (GCurrent[locBreak][i] < 1 - 1e-6)&(GCurrent[locBreak][i] > 1e-6)
                        θgainTemp += GCurrent[locBreak][i]*disData[locBreak].d[i];
                    end
                else
                    if (GCurrent[locBreak][i] < 1 - 1e-6)&(GCurrent[locBreak][i] > 1e-6)
                        θgainTemp += min(1 - GCurrent[locBreak][i],GCurrent[locBreak][i])*disData[locBreak].d[i];
                    end
                end
            end
            if θgainTemp > θgain
                θgain = θgainTemp;
                lGFracInd = i;
            end
        end
        # initialize the breakPoints dictionary
        if lGFracInd != -1
            #locBreak = Int64(floor((GFrac[lGFracInd][1]*fracBreak + GFrac[lGFracInd][2]*(1 - fracBreak))));
            divSet1,divDet1,divSet2,divDet2 = breakDiv(pData,disData,H,divSet,divDet,lGFracInd,locBreakH,distanceDict);
            divSet1,divDet1 = divExploit(pData,disData,H,divSet1,divDet1,distanceDict);
            divSet1,divDet1 = splitPrepld2(pData,disData,Ω,H,GCurrent,tCurrent,divSet1,divDet1,θCurrent,θIntCurrent,nSplit);
            #divSet1,divDet1 = splitPrepSmart2(pData,disData,Ω,H,GCurrent,tCurrent,divSet1,divDet1,θCurrent,θIntCurrent,nSplit)

            divSet2,divDet2 = divExploit(pData,disData,H,divSet2,divDet2,distanceDict);
            divSet2,divDet2 = splitPrepld2(pData,disData,Ω,H,GCurrent,tCurrent,divSet2,divDet2,θCurrent,θIntCurrent,nSplit);
            #divSet2,divDet2 = splitPrepSmart2(pData,disData,Ω,H,GCurrent,tCurrent,divSet2,divDet2,θCurrent,θIntCurrent,nSplit)
            returnSet = [[divSet1,divDet1],[divSet2,divDet2]];
            # select the cuts
            if cutSelOpt
                cutSel = examineCuts_count_3(pData,disData,Ω,cutSetNew,divSet,tCurrent,xCurrent,θCurrent,yCurrent,IJPair,IPPair);
                cutSetRe = selectCuts3(cutSetNew,cutSel);
            else
                cutSetRe = deepcopy(cutSetNew);
            end
        else
            # if all i's having binary G's, we reach optimum for this node, ub = lb
            returnSet = [];
            cutSetRe = [];
        end
    elseif (mpStatus == :UserLimit)
        returnNo = getobjectivebound(mp);
        if !(isnan(mpObj))
            for i in pData.II
                tCurrent[i] = getvalue(mp[:t][i]);
                for j in pData.Ji[i]
                    xCurrent[i,j] = getvalue(mp[:x][i,j]);
                end
                for par in 1:length(divSet[i])
                    yCurrent[i,par] = getvalue(mp[:y][i,par]);
                end
            end
            subInfo = pmap(wp,ω -> sub_divT(pData,disData[ω],ω,tCurrent,xCurrent,yCurrent,divSet,H,lDict,2),Ω);
            GCurrent = [subInfo[ω][2] for ω in Ω];
            θCurrent = [subInfo[ω][1] for ω in Ω];
            ubCurrent,θIntCurrent = ubCalP(pData,disData,Ω,xCurrent,tCurrent,Tmax1,1,wp);
            # branch
            θDiff = [θIntCurrent[ω] - θCurrent[ω] for ω in Ω];
            θDiffPerm = sortperm(θDiff,rev = true);
            locBreak = θDiffPerm[1];
            locBreakH = HRev[disData[locBreak].H];
            θgain = 0;
            lGFracInd = -1;
            for i in pData.II
                θgainTemp = 0;
                if i != 0
                    if tCurrent[i] > H[locBreakH] + 1e-6
                        if (GCurrent[locBreak][i] < 1 - 1e-6)&(GCurrent[locBreak][i] > 1e-6)
                            θgainTemp += (1 - GCurrent[locBreak][i])*disData[locBreak].d[i];
                        end
                    elseif tCurrent[i] < H[locBreakH] - 1e-6
                        if (GCurrent[locBreak][i] < 1 - 1e-6)&(GCurrent[locBreak][i] > 1e-6)
                            θgainTemp += GCurrent[locBreak][i]*disData[locBreak].d[i];
                        end
                    else
                        if (GCurrent[locBreak][i] < 1 - 1e-6)&(GCurrent[locBreak][i] > 1e-6)
                            θgainTemp += min(1 - GCurrent[locBreak][i],GCurrent[locBreak][i])*disData[locBreak].d[i];
                        end
                    end
                end
                if θgainTemp > θgain
                    θgain = θgainTemp;
                    lGFracInd = i;
                end
            end
            # initialize the breakPoints dictionary
            if lGFracInd != -1
                #locBreak = Int64(floor((GFrac[lGFracInd][1]*fracBreak + GFrac[lGFracInd][2]*(1 - fracBreak))));
                divSet1,divDet1,divSet2,divDet2 = breakDiv(pData,disData,H,divSet,divDet,lGFracInd,locBreakH,distanceDict);
                divSet1,divDet1 = divExploit(pData,disData,H,divSet1,divDet1,distanceDict);
                divSet1,divDet1 = splitPrepld2(pData,disData,Ω,H,GCurrent,tCurrent,divSet1,divDet1,θCurrent,θIntCurrent,nSplit);
                #divSet1,divDet1 = splitPrepSmart2(pData,disData,Ω,H,GCurrent,tCurrent,divSet1,divDet1,θCurrent,θIntCurrent,nSplit)

                divSet2,divDet2 = divExploit(pData,disData,H,divSet2,divDet2,distanceDict);
                divSet2,divDet2 = splitPrepld2(pData,disData,Ω,H,GCurrent,tCurrent,divSet2,divDet2,θCurrent,θIntCurrent,nSplit);
                #divSet2,divDet2 = splitPrepSmart2(pData,disData,Ω,H,GCurrent,tCurrent,divSet2,divDet2,θCurrent,θIntCurrent,nSplit)
                returnSet = [[divSet1,divDet1],[divSet2,divDet2]];
                if cutSelOpt
                    cutSel = examineCuts_count_3(pData,disData,Ω,cutSetNew,divSet,tCurrent,xCurrent,θCurrent,yCurrent,IJPair,IPPair);
                    cutSetRe = selectCuts3(cutSetNew,cutSel);
                else
                    cutSetRe = deepcopy(cutSetNew);
                end
            else
                # if all i's having binary G's, we reach optimum for this node, ub = lb
                returnSet = [];
                cutSetRe = [];
            end
        else
            # branch only + equal branching rules
            # select the activity with the largest variance in the past solution to branch
            lGFracInd = -1;
            maxVar = 0;
            for i in pData.II
                tHistory = [tcoreList[ll][i] for ll in 1:length(tcoreList)];
                if var(tHistory) > maxVar
                    lGFracInd = i;
                    maxVar = var(tHistory);
                end
            end
            earliest0 = 0;
            latest0 = length(H) - 1;
            for par in 1:length(divSet[lGFracInd])
                if divDet[lGFracInd][par] == 1
                    earliest0 = divSet[lGFracInd][par].endH;
                end
                if (divDet[lGFracInd][par] == -1)&(divSet[lGFracInd][par].startH < latest0)
                    latest0 = divSet[lGFracInd][par].startH;
                end
            end
            if latest0 - earliest0 > 1
                locBreak = Int64(floor((earliest0 + latest0)/2));
                locBreakH = HRev[disData[locBreak].H];
                divSet1,divDet1,divSet2,divDet2 = breakDiv(pData,disData,H,divSet,divDet,lGFracInd,locBreakH,distanceDict);
                divSet1,divDet1 = divExploit(pData,disData,H,divSet1,divDet1,distanceDict);
                divSet2,divDet2 = divExploit(pData,disData,H,divSet2,divDet2,distanceDict);
                returnSet = [[divSet1,divDet1],[divSet2,divDet2]];
                cutSetRe = [];
            else
                returnSet = [[divSet,divDet]];
                cutSetRe = [];
            end
        end
    else
        returnNo = Inf;
        returnSet = [];
        cutSetRe = [];
    end
    return returnNo,cutSetRe,returnSet,tbest,xbest,minimum(ubCostList),tcoreNew,xcoreNew,ubcoreNew;
end


function runPara_Share(treeList,cutList,tcoreList,xcoreList,ubcoreList,ubCost,tbest,xbest,batchNo,noTh,ϵ = 1e-2,nSplit = 5,noPa = 1,cutSelOpt = true)
    # separate the workers to main processors and workers
    npList = workers()[1:batchNo];
    global noMo = div(noThreads,batchNo);
    # noPa = noMo - noTh;
    wpDict = Dict();
    for npi in 1:length(npList)
        wpDict[npList[npi]] = workers()[(batchNo + (npi - 1)*noPa + 1):(batchNo + npi*noPa)];
    end
    lbOverAll = -Inf;
    timeDict = Dict();
    lbDict = Dict();
    lbDict[1] = -Inf;

    @sync begin
        for ip in 1:length(npList)
            p = npList[ip];
            @async begin
                while true
                    # if all nodes are processed and no nodes are being processed, exit
                    boolFinished = true;
                    if (ubCost - lbOverAll)/ubCost >= ϵ
                        for l in 1:length(treeList)
                            if treeList[l][3] != 1
                                boolFinished = false;
                            end
                        end
                    end
                    println("-------------",boolFinished," ",[treeList[l][3] for l in 1:length(treeList)],[treeList[l][1] for l in 1:length(treeList)],"-------------");
                    if boolFinished
                        println("break");
                        break
                    end
                    openNodes = [(treeList[l][1],l) for l in 1:length(treeList) if treeList[l][3] == -1];
                    if openNodes != []
                        selectNode = sort(openNodes, by = x -> x[1])[1][2];
                        if (treeList[selectNode][1] < ubCost) & ((ubCost - lbDict[selectNode])/ubCost >= ϵ)
                            println("On core: ",p," processing node: ",selectNode," lower bound is: ",treeList[selectNode][1]," upper bound is: ",minimum(ubcoreList));
                            treeList[selectNode][3] = 0;
                            cutData = cutList[treeList[selectNode][2]];
                            divData = [treeList[id][4] for id in treeList[selectNode][2]];
                            tic();
                            # mpSolveInfo = remotecall_fetch(solveMP_para_Share,p,[cutData,divData,treeList[selectNode][4],tcoreList,xcoreList,ubcoreList,ubCost,
                            #     tbest,xbest,noTh,wpDict[p],nSplit,pData,disData,lDict,H,allSucc,distanceDict]);
                            mpSolveInfo = remotecall_fetch(solveMP_para_Share,p,[cutData,divData,treeList[selectNode][4],tcoreList,xcoreList,ubcoreList,ubCost,
                                tbest,xbest,noTh,wpDict[p],nSplit,treeList[selectNode][5],cutSelOpt]);
                            timeDict[selectNode] = toc();
                            # update the cutList with the added cuts and two new nodes
                            # update the cutSet
                            # return returnNo,cutSet,returnSet,tbest,xbest,minimum(ubCostList)
                            treeList[selectNode][3] = 1;
                            # the problem is solved and a lower bound is generated
                            lbDict[selectNode] = mpSolveInfo[1];
                            if mpSolveInfo[1] < Inf
                                append!(tcoreList,mpSolveInfo[7]);
                                append!(xcoreList,mpSolveInfo[8]);
                                append!(ubcoreList,mpSolveInfo[9]);
                                if mpSolveInfo[6] < ubCost
                                    ubCost = mpSolveInfo[6];
                                    tbest = mpSolveInfo[4];
                                    xbest = mpSolveInfo[5];
                                end
                                # update the current node cut
                                cutList[selectNode] = mpSolveInfo[2];
                                # push the branched nodes
                                if (mpSolveInfo[1] < ubCost)
                                    ancestorTemp = deepcopy(treeList[selectNode][2]);
                                    push!(ancestorTemp,selectNode);
                                    rlTemp = treeList[selectNode][5];
                                    if mpSolveInfo[1] > maximum([treeList[l][1] for l in ancestorTemp])
                                        lbNode = mpSolveInfo[1];
                                    else
                                        lbNode = maximum([treeList[l][1] for l in ancestorTemp]);
                                        rlTemp = rlTemp * 2;
                                    end
                                    for newN in 1:length(mpSolveInfo[3])
                                        push!(treeList,[lbNode,ancestorTemp,-1,mpSolveInfo[3][newN],rlTemp]);
                                        push!(cutList,[]);
                                        # before solving the node, the node's lower bound should be its father node's lb value
                                        lbDict[length(treeList)] = lbNode;
                                    end
                                end
                            end
                            if [treeList[l][1] for l in 1:length(treeList) if treeList[l][3] != 1] != []
                                lbOverAll = minimum([treeList[l][1] for l in 1:length(treeList) if treeList[l][3] != 1]);
                            end
                        end
                        treeList[selectNode][3] = 1;
                    else
                        println("**************Worker $(p) waiting for open nodes.**************");
                        remotecall_fetch(sleep, p, 10);
                    end
                end
            end
        end
    end

    # search the tree for all leaf nodes to obtain the lower bound
    treeStruct = Dict();
    for l in 1:length(treeList)
        treeStruct[l] = [];
        for item in treeList[l][2]
            push!(treeStruct[item],l);
        end
    end
    minLB = Inf;
    for l in 1:length(treeList)
        if treeStruct[l] == []
            if l in keys(lbDict)
                if lbDict[l] < minLB
                    minLB = lbDict[l];
                end
            end
        end
    end

    return tbest,xbest,ubCost,minLB,timeDict,treeList;
end

function runPara_Series_Share(treeList,cutList,tcoreList,xcoreList,ubcoreList,ubCost,tbest,xbest,noTh,nSplit = 5)
    # separate the workers to main processors and workers
    boolFinished = true;
    while boolFinished
        # if all nodes are processed and no nodes are being processed, exit
        println("-------------",boolFinished," ",[treeList[l][3] for l in 1:length(treeList)],[treeList[l][1] for l in 1:length(treeList)],"-------------");
        openNodes = [(treeList[l][1],l) for l in 1:length(treeList) if treeList[l][3] == -1];
        if openNodes != []
            selectNode = sort(openNodes, by = x -> x[1])[1][2];
            if treeList[selectNode][1] < ubCost
                println("Processing node: ",selectNode," lower bound is: ",treeList[selectNode][1]," upper bound is: ",minimum(ubcoreList));
                treeList[selectNode][3] = 0;
                cutData = cutList[treeList[selectNode][2]];
                divData = [treeList[id][4] for id in treeList[selectNode][2]];
                mpSolveInfo = solveMP_para_Share([cutData,divData,treeList[selectNode][4],tcoreList,xcoreList,ubcoreList,ubCost,tbest,xbest,noTh,workers(),nSplit]);
                # update the cutList with the added cuts and two new nodes
                # update the cutSet
                # return returnNo,cutSet,returnSet,tbest,xbest,minimum(ubCostList)
                treeList[selectNode][3] = 1;
                if mpSolveInfo[1] > -Inf
                    append!(tcoreList,mpSolveInfo[7]);
                    append!(xcoreList,mpSolveInfo[8]);
                    append!(ubcoreList,mpSolveInfo[9]);
                    if mpSolveInfo[6] < ubCost
                        ubCost = mpSolveInfo[6];
                        tbest = mpSolveInfo[4];
                        xbest = mpSolveInfo[5];
                    end
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
                if [treeList[l][1] for l in 1:length(treeList) if treeList[l][3] != 1] != []
                    lbOverAll = minimum([treeList[l][1] for l in 1:length(treeList) if treeList[l][3] != 1]);
                    if (ubCost - lbOverAll)/ubCost < ϵ
                        boolFinished = false;
                    end
                end
            end
        else
            boolFinished = false;
        end
    end

    return tbest,xbest,ubCost,lbOverAll;
end

function partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,batchNo,noTh,noPa,ϵ = 1e-2,nSplit = 5,roundLimit = 1000,cutSelOpt = true)
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
    DShare = SharedArray{Float64,1}(length(pData.II));
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
    # data141 = load("14_test1_ubData.jld");
    # ubextList,tHList,ubInc,tbest,xbest,θbest,textList,xextList = data141["data"];
    ubextList,tHList,ubInc,tbest,xbest,θbest,textList,xextList = iniPart(pData,disData,Ω,sN,MM,1,noThreads);
    lbCost = -Inf;
    lbCostList = [];
    global ubCost = ubInc;
    tcoreList = deepcopy(textList);
    xcoreList = deepcopy(xextList);
    ubcoreList = deepcopy(ubextList);

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

    # divInfoShare = convertDiv(divSet,divDet);

    # set up a tree list
    global treeList = [];
    global cutList = [];
    push!(treeList,[lbCost,[],-1,[divSet,divDet],roundLimit]); # the empty set is the list of predecessors of the current node
    push!(cutList,[]);

    global lbOverAll = -Inf;
    # transfer the data back to everywhere
    tic();
    tbest,xbest,ubCost,lbOverAll,timeIter,treeList = runPara_Share(treeList,cutList,textList,xextList,ubextList,ubCost,tbest,xbest,batchNo,noTh,ϵ,nSplit,noPa,cutSelOpt);
    decompTime = toc();

    # need a cut selection process within the callback
    return tbest,xbest,ubCost,lbOverAll,timeIter,treeList,decompTime;
end
