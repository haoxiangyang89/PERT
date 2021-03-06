function subPara1(pData,disData,Ω,tbest,xbest,ybest,divSet,H,lDict,wp = CachingPool(workers()))
    θList = pmap(ω -> sub_divT(pData,disData[ω],ω,tbest,xbest,ybest,divSet,H,lDict),wp,Ω);
    return θList;
end

function subPara(pData,disData,Ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore,wp = CachingPool(workers()))
    dataList = pmap(ω -> sub_divTDualT2(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore), wp, Ω);
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

function solveMP_para(data)
    # input: [cutData,cutCurrent,tcoreList,xcoreList,ubcoreList,ubCost,tbest,xbest,noTh,wpList]
    pData = data[13];
    disData = data[14];
    divSet,divDet = data[3];
    divData = data[2];
    cutSet = data[1];           # historical cuts
    IJPair = [(i,j) for i in pData.II for j in pData.Ji[i]];
    IPPair = [(i,par) for i in pData.II for par in 1:length(divSet[i])];
    tcoreList = data[4];
    xcoreList = data[5];
    ubcoreList = data[6];
    ubCost = data[7];
    tbest = data[8];
    xbest = data[9];
    noTh = data[10];
    wp = CachingPool(data[11]);
    nSplit = data[12];
    lDict = data[15];
    H = data[16];
    allSucc = data[17];
    distanceDict = data[18];
    Ω = 1:length(disData);

    Tmax1 =lDict[0];
    GList = [];
    tcoreNew = [];
    xcoreNew = [];
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
            #that = SharedArray{Float64,1}((length(pData.II)));
            that = Dict();
            #xhat = SharedArray{Float64,1}((length(IJPair)));
            xhat = Dict();
            #θhat = SharedArray{Float64,1}((length(Ω)));
            #yhat = SharedArray{Float64,1}((length(IPPair)));
            yhat = Dict();
            θhat = Dict();
            # obtain the solution at the current node
            for i in pData.II
                #that[findfirst(x -> x==i, pData.II)] = getvalue(t[i]);
                that[i] = getvalue(t[i]);
                for j in pData.Ji[i]
                    #xhat[findfirst(x -> x == (i,j), IJPair)] = getvalue(x[i,j]);
                    xhat[i,j] = getvalue(x[i,j]);
                end
                for par in 1:length(divSet[i])
                    #yhat[findfirst(x -> x == (i,par), IPPair)] = round(getvalue(y[i,par]));
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
            for ω in Ω
                infeasBool = false;
                if length(dataList[ω]) == 3
                    push!(tError,that);
                    push!(xError,xhat);
                    push!(yError,yhat);
                    push!(ErrorωList,ω);
                    push!(tcoreError,tcore);
                    push!(xcoreError,xcore);
                    push!(ycoreError,ycore);
                    infeasBool = true;
                elseif length(dataList[ω]) == 6
                    push!(tUnbounded,that);
                    push!(xUnbounded,xhat);
                    push!(yUnbounded,yhat);
                    push!(UnboundedωList,ω);
                    push!(tcoreUnbounded,tcore);
                    push!(xcoreUnbounded,xcore);
                    push!(ycoreUnbounded,ycore);
                end
                if infeasBool
                    return JuMP.StopTheSolver;
                end
            end
            cutScen = [ω for ω in Ω if dataList[ω][4] - θhat[findfirst(x -> x==ω,Ω)] > 1e-4*θhat[findfirst(x -> x==ω,Ω)]];
            # πSet = SharedArray{Float64,2}((length(pData.II),length(cutScen)));
            # λSet = SharedArray{Float64,2}((length(IJPair),length(cutScen)));
            # γSet = SharedArray{Float64,2}((length(IPPair),length(cutScen)));
            # vSet = SharedArray{Float64,1}((length(cutScen)));
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
                        πSet[findfirst(x -> x==i, pData.II),ωi] = dataList[ω][1][i];
                    else
                        πSet[findfirst(x -> x==i, pData.II),ωi] = 0;
                        if dataList[ω][1][i] < 0
                            vSet[ωi] += dataList[ω][1][i];
                        end
                    end
                    for j in pData.Ji[i]
                        vSet[ωi] -= dataList[ω][2][i,j]*xhat[i,j];
                        if abs(dataList[ω][2][i,j]) >= 1e-5
                            λSet[findfirst(x -> x == (i,j), IJPair),ωi] = dataList[ω][2][i,j];
                        else
                            λSet[findfirst(x -> x == (i,j), IJPair),ωi] = 0;
                            if dataList[ω][2][i,j] < 0
                                vSet[ωi] += dataList[ω][2][i,j];
                            end
                        end
                    end
                    for par in 1:length(divSet[i])
                        vSet[ωi] -= dataList[ω][3][i,par]*yhat[i,par];
                        if abs(dataList[ω][3][i,par]) >= 1e-5
                            γSet[findfirst(x -> x == (i,par), IPPair),ωi] = dataList[ω][3][i,par];
                        else
                            γSet[findfirst(x -> x == (i,par), IPPair),ωi] = 0;
                            if dataList[ω][3][i,par] < 0
                                vSet[ωi] += dataList[ω][3][i,par];
                            end
                        end
                    end
                end
                @lazyconstraint(cb, θ[ω] >= vSet[ωi] + sum(πSet[findfirst(x -> x==i, pData.II),ωi]*t[i] for i in pData.II) +
                    sum(sum(λSet[findfirst(x -> x == (i,j), IJPair),ωi]*x[i,j] for j in pData.Ji[i]) for i in pData.II) +
                    sum(sum(γSet[findfirst(x -> x == (i,par), IPPair),ωi]*y[i,par] for par in 1:length(divSet[i])) for i in pData.II));
            end
            newCuts = [cutScen,πSet,λSet,γSet,vSet];
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
    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-8, FeasibilityTol = 1e-8, Method = 1, Threads = noTh, Cutoff = ubCost));
    #mp = Model(solver = GurobiSolver(IntFeasTol = 1e-8, FeasibilityTol = 1e-8, Method = 1, NumericFocus = 3, Threads = noTh, Cutoff = ubCost));
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
                @constraint(mp, θ[ω] >= vk + sum(πk[findfirst(x -> x==i, pData.II)]*mp[:t][i] +
                    sum(λk[findfirst(x -> x == (i,j), IJPair)]*mp[:x][i,j] for j in pData.Ji[i]) +
                    sum(γk[findfirst(x -> x == (i,par), IPPairPrev)]*(sum(mp[:y][i,parNew] for parNew in revDict[i][par]))
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
        for i in pData.II
            tCurrent[i] = getvalue(mp[:t][i]);
            for j in pData.Ji[i]
                xCurrent[i,j] = getvalue(mp[:x][i,j]);
            end
            for par in 1:length(divSet[i])
                yCurrent[i,par] = getvalue(mp[:y][i,par]);
            end
        end
        # for ω in Ω
        #     θCurrent[ω] = getvalue(mp[:θ][ω]);
        # end
        θCurrent = pmap(ω -> sub_divT(pData,disData[ω],ω,tCurrent,xCurrent,yCurrent,divSet,H,lDict),wp,Ω);
        ubCurrent,θIntCurrent = ubCalP(pData,disData,Ω,xCurrent,tCurrent,Tmax1,1,wp);
        # branch
        GCurrent = GList[length(GList)];
        θDiff = [θIntCurrent[ω] - θCurrent[ω] for ω in Ω];
        θDiffPerm = sortperm(θDiff,rev = true);
        locBreak = θDiffPerm[1];
        θgain = 0;
        lGFracInd = -1;
        for i in pData.II
            θgainTemp = 0;
            if i != 0
                if tCurrent[i] > H[locBreak] + 1e-6
                    if (GCurrent[locBreak][i] < 1 - 1e-6)&(GCurrent[locBreak][i] > 1e-6)
                        θgainTemp += (1 - GCurrent[locBreak][i])*disData[locBreak].d[i];
                    end
                elseif tCurrent[i] < H[locBreak] - 1e-6
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
            divSet1,divDet1,divSet2,divDet2 = breakDiv(pData,disData,H,divSet,divDet,lGFracInd,locBreak,distanceDict);
            divSet1,divDet1 = divExploit(pData,disData,H,divSet1,divDet1,distanceDict);
            divDict = Dict();
            for i in pData.II
                divDict[i] = [];
                θDiv = [];
                for ω in Ω
                    currentpar = -1;
                    for par in 1:length(divSet1[i])
                        if (ω > divSet1[i][par].startH) & (ω < divSet1[i][par].endH)
                            currentpar = par;
                        end
                    end
                    if currentpar != -1
                        if divDet1[i][currentpar] == 0
                            if i != 0
                                if tCurrent[i] > H[ω] + 1e-6
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,((1 - GCurrent[ω][i])*disData[ω].d[i],ω));
                                    end
                                elseif tCurrent[i] < H[ω] - 1e-6
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,(GCurrent[ω][i]*disData[ω].d[i],ω));
                                    end
                                else
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,(min(1 - GCurrent[ω][i],GCurrent[ω][i])*disData[ω].d[i],ω));
                                    end
                                end
                            else
                                if tCurrent[i] > H[ω] + 1e-6
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,(1 - GCurrent[ω][i],ω));
                                    end
                                elseif tCurrent[i] < H[ω] - 1e-6
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,(GCurrent[ω][i],ω));
                                    end
                                else
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,(min(1 - GCurrent[ω][i],GCurrent[ω][i]),ω));
                                    end
                                end
                            end
                        end
                    end
                end
                θDivperm = sortperm(θDiv,rev = true);
                if θDiv != []
                    # if there are fractional solutions
                    for n in 1:(min(nSplit,length(θDiv)))
                        ωn = θDiv[θDivperm[n]][2];
                        push!(divDict[i],ωn);
                    end
                end
            end
            divSet1,divDet1 = splitAny(divSet1,divDet1,Ω,divDict);

            divSet2,divDet2 = divExploit(pData,disData,H,divSet2,divDet2,distanceDict);
            divDict = Dict();
            for i in pData.II
                divDict[i] = [];
                θDiv = [];
                for ω in Ω
                    currentpar = -1;
                    for par in 1:length(divSet2[i])
                        if (ω > divSet2[i][par].startH) & (ω < divSet2[i][par].endH)
                            currentpar = par;
                        end
                    end
                    if currentpar != -1
                        if divDet2[i][currentpar] == 0
                            if i != 0
                                if tCurrent[i] > H[ω] + 1e-6
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,((1 - GCurrent[ω][i])*disData[ω].d[i],ω));
                                    end
                                elseif tCurrent[i] < H[ω] - 1e-6
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,(GCurrent[ω][i]*disData[ω].d[i],ω));
                                    end
                                else
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,(min(1 - GCurrent[ω][i],GCurrent[ω][i])*disData[ω].d[i],ω));
                                    end
                                end
                            else
                                if tCurrent[i] > H[ω] + 1e-6
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,(1 - GCurrent[ω][i],ω));
                                    end
                                elseif tCurrent[i] < H[ω] - 1e-6
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,(GCurrent[ω][i],ω));
                                    end
                                else
                                    if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)
                                        push!(θDiv,(min(1 - GCurrent[ω][i],GCurrent[ω][i]),ω));
                                    end
                                end
                            end
                        end
                    end
                end
                θDivperm = sortperm(θDiv,rev = true);
                if θDiv != []
                    # if there are fractional solutions
                    for n in 1:(min(nSplit,length(θDiv)))
                        ωn = θDiv[θDivperm[n]][2];
                        push!(divDict[i],ωn);
                    end
                end
            end
            divSet2,divDet2 = splitAny(divSet2,divDet2,Ω,divDict);
            returnSet = [[divSet1,divDet1],[divSet2,divDet2]];
        else
            # if all i's having binary G's, we reach optimum for this node, ub = lb
            returnSet = [];
        end
    else
        returnNo = -Inf;
        returnSet = [];
    end
    return returnNo,cutSetNew,returnSet,tbest,xbest,minimum(ubCostList),tcoreNew,xcoreNew,ubcoreNew;
end

function partSolve_BB_para(pData,disData,Ω,sN,MM,noThreads,batchNo = 5,ϵ = 1e-2,nSplit = 5)
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
    #ubextList,tHList,ubInc,tbest,xbest,θbest,textList,xextList = iniPart(pData,disData,Ω,sN,MM,1,noThreads);
    data141 = load("14_test1_ubData.jld");
    ubextList,tHList,ubInc,tbest,xbest,θbest,textList,xextList = data141["data"];
    lbCost = -Inf;
    lbCostList = [];
    global ubCost = ubInc;

    tcoreList = deepcopy(textList);
    xcoreList = deepcopy(xextList);
    ubcoreList = deepcopy(ubextList);
    ycoreList = [];
    errorList = [];
    GList = [];

    brInfo = precludeRelNew(pData,H,ubCost);

    # initialize cutSet and divSet
    cutSet = [];
    divSet = Dict();
    divDet = Dict();
    for i in pData.II
        set1 = [h for h in HΩ if brInfo[findfirst(x -> x==i, pData.II),h] == 1];
        setn1 = [h for h in HΩ if brInfo[findfirst(x -> x==i, pData.II),h] == -1];

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

    # set up a tree list
    global treeList = [];
    global cutList = [];
    push!(treeList,[lbCost,[],-1,[divSet,divDet]]); # the empty set is the list of predecessors of the current node
    push!(cutList,[]);

    global lbOverAll = -Inf;

    # separate the workers to main processors and workers
    npList = workers()[1:batchNo];
    global noTh = div(noThreads,batchNo);
    noPa = noTh - 1;
    wpDict = Dict();
    for npi in 1:length(npList)
        wpDict[npList[npi]] = workers()[(batchNo + (npi - 1)*noPa + 1):(batchNo + npi*noPa)];
    end
    keepIter = true;
    lbOverAll = -Inf;
    timeDict = Dict();

    @sync begin
        for ip in 1:length(npList)
            p = npList[ip];
            tError = @async begin
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
                        if treeList[selectNode][1] < ubCost
                            println("On core: ",p," processing node: ",selectNode," lower bound is: ",treeList[selectNode][1]," upper bound is: ",minimum(ubcoreList));
                            treeList[selectNode][3] = 0;
                            cutData = cutList[treeList[selectNode][2]];
                            divData = [treeList[id][4] for id in treeList[selectNode][2]];
                            tempTimer = time();
                            mpSolveInfo = remotecall_fetch(solveMP_para,p,[cutData,divData,treeList[selectNode][4],tcoreList,xcoreList,ubcoreList,ubCost,
                                tbest,xbest,noTh,wpDict[p],nSplit,pData,disData,lDict,H,allSucc,distanceDict]);
                            timeDict[selectNode] = time() - tempTimer;
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
            wait(tError);
        end
    end

    # need a cut selection process within the callback
    return tbest,xbest,ubCost,lbOverAll,timeDict;
end
