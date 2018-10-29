function partitionSolve(pData,disData,ϵ = 0.01,tightenBool = 0, cutThreshold = 10)
    # process to solve the PERT problem
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
    ubCost = ubdet;
    lbCost = -Inf;
    lbCostList = [];
    ubCostList = [];
    # set up the counter for being tight
    cutSel = Dict();
    cutyn = [];
    cutynRec = [];

    while (ubCost - lbCost)/ubCost > ϵ
        keepIter = true;
        tlb = Dict();
        xlb = Dict();
        θlb = Dict();
        ylb = Dict();
        dataList = [];
        mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax,0,cutyn);
        # if perform the bound tightening process
        if tightenBool == 1
            mpTemp = copy(mp);
            ubInfo,lbInfo = obtainBds(pData,disData,Ω,mpTemp,ubCost);
            mp = updateMaster(mp,ubInfo,lbInfo);
            divSet,divDet = revisePar(pData,disData,divSet,divDet,ubInfo,lbInfo);
            mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax);
        elseif tightenBool == 2
            mpTemp = copy(mp);
            divDet = obtainDet(pData,disData,Ω,mpTemp,ubCost,divSet,divDet);
            mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax);
        end
        while keepIter
            mp.solver = GurobiSolver();
            solve(mp);
            # obtain the solution
            that = Dict();
            xhat = Dict();
            θhat = Dict();
            yhat = Dict();
            for i in pData.II
                that[i] = getvalue(mp[:t][i]);
                for j in pData.Ji[i]
                    xhat[i,j] = getvalue(mp[:x][i,j]);
                end
                for par in 1:length(divSet[i])
                    yhat[i,par] = getvalue(mp[:y][i,par]);
                end
            end
            for ω in Ω
                θhat[ω] = getvalue(mp[:θ][ω]);
            end
            lbCost = getobjectivevalue(mp);
            if lbCostList != []
                if lbCost < maximum(lbCostList)
                    cutThreshold += 5;
                end
            end
            push!(lbCostList,lbCost);
            # examine how many cuts are tight at this solution, update the cutSel
            cutSel,cutyn = examineCuts_count(disData,Ω,cutSel,cutSet,that,xhat,θhat,yhat,cutThreshold);
            push!(cutynRec,length(cutyn));
            # update the master with new cutyn
            mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax,0,cutyn);

            # generate cuts
            lbPrev = lbCost;
            πdict = Dict();
            λdict = Dict();
            γdict = Dict();
            vk = Dict();
            θInt = Dict();
            ubTemp,θInt = ubCalP(pData,disData,Ω,xhat,that,Tmax1,1);
            if ubCost > ubTemp
                ubCost = ubTemp;
                tbest = copy(that);
                xbest = copy(xhat);
            end
            push!(ubCostList,ubCost);
            dataList = pmap(ω -> sub_divT(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict), Ω);
            for ω in Ω
                πdict[ω] = dataList[ω][1];
                λdict[ω] = dataList[ω][2];
                γdict[ω] = dataList[ω][3];
                vk[ω] = dataList[ω][4];
            end
            ωTightCounter = 0;
            cutDual = [];
            for ω in Ω
                if vk[ω] - θhat[ω] > 1e-4*θhat[ω]
                    push!(cutDual,[ω,vk[ω],πdict[ω],λdict[ω],γdict[ω]]);
                    mp = addtxyCut(pData,ω,mp,πdict[ω],λdict[ω],γdict[ω],vk[ω],that,xhat,yhat,divSet);
                else
                    ωTightCounter += 1;
                end
            end
            if ωTightCounter == length(Ω)
                keepIter = false;
                for i in pData.II
                    tlb[i] = that[i];
                    for j in pData.Ji[i]
                        xlb[i,j] = xhat[i,j];
                    end
                    for par in 1:length(divSet[i])
                        ylb[i,par] = yhat[i,par];
                    end
                end
                for ω in Ω
                    θlb[ω] = θhat[ω];
                end
            else
                push!(cutSet,[[that,xhat,yhat,divSet],cutDual]);
                for ω in Ω
                    cutSel[length(cutSet),ω] = 0;
                end
            end
        end

        # need to come up with a rule to partition: gradient descent like binary search
        # check θInt vs. θhat: why the lower bound and the upper bound do not converge quickly --->
        # use the sub problem solution G to learn the b&b
        # also need to think up a way to tightening the cuts for each partition
        GFrac = Dict();
        for i in pData.II
            GFraciList = [ω for ω in Ω if (dataList[ω][5][i] < 1 - 1e-6)&(dataList[ω][5][i] > 1e-6)];
            if GFraciList != []
                GFrac[i] = [minimum(GFraciList),maximum(GFraciList)];
            else
                GFrac[i] = [];
            end
        end
        # create new partition
        newPartition = [];
        for i in pData.II
            if GFrac[i] != []
                newItem = (i,Int(floor((GFrac[i][1] + GFrac[i][2])/2)));
                push!(newPartition,newItem);
            end
        end
        divSet,divDet = splitPar(divSet,divDet,newPartition);
    end

    return tbest,xbest,lbCost,ubCost;
end

function limYselection(pData,H,tcurr,divSet,radius)
    # identify the partition within radius to the current t
    yLim = Dict();
    for i in pData.II
        # if tcurr[i] in divSet[par]
        yLim[i] = [];
        for par in 1:length(divSet[i])
            if (tcurr[i] <= H[divSet[i][par].endH])&(tcurr[i] >= H[divSet[i][par].startH])
                # for every partition in the radius
                if radius != Inf
                    startPar = max(par - radius,1);
                    endPar = min(par + radius,length(divSet[i]));
                else
                    startPar = 1;
                    endPar = length(divSet[i]);
                end
                for par1 in startPar:endPar
                    if !(par1 in yLim[i])
                        push!(yLim[i],par1);
                    end
                end
            end
        end
    end
    return yLim;
end

function partitionSolve_yLim(pData,disData,distanceDict,allSucc,ϵ = 0.01,tightenBool = 0, cutThreshold = 10,radius = 1,partOption = 2)
    # process to solve the PERT problem
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
    ubCost = ubdet;
    lbCost = -Inf;
    lbCostList = [];
    ubCostList = [];
    # set up the counter for being tight
    cutSel = Dict();
    cutyn = [];
    cutynRec = [];
    yLim = limYselection(pData,H,tdet,divSet,radius);
    tcoreList = [];
    xcoreList = [];
    ycoreList = [];

    while (ubCost - lbCost)/ubCost > ϵ
        keepIter = true;
        tlb = Dict();
        xlb = Dict();
        θlb = Dict();
        ylb = Dict();
        dataList = [];
        mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax,distanceDict,allSucc,1,yLim,0,cutyn);
        # if perform the bound tightening process
        if tightenBool == 1
            mpTemp = copy(mp);
            ubInfo,lbInfo = obtainBds(pData,disData,Ω,mpTemp,ubCost);
            mp = updateMaster(mp,ubInfo,lbInfo);
            divSet,divDet = revisePar(pData,disData,divSet,divDet,ubInfo,lbInfo);
            mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax,distanceDict,allSucc,1,yLim);
        elseif tightenBool == 2
            mpTemp = copy(mp);
            divDet = obtainDet(pData,disData,Ω,mpTemp,ubCost,divSet,divDet);
            mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax,distanceDict,allSucc,1,yLim);
        end

        tRec = Dict();
        for i in pData.II
            tRec[i] = [];
        end
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
        tic();
        while keepIter
            mp.solver = CplexSolver(CPX_PARAM_EPINT = 1e-9,CPX_PARAM_EPRHS = 1e-9);
            solve(mp);
            # obtain the solution
            that = Dict();
            xhat = Dict();
            θhat = Dict();
            yhat = Dict();
            for i in pData.II
                that[i] = getvalue(mp[:t][i]);
                push!(tRec[i],that[i]);
                for j in pData.Ji[i]
                    xhat[i,j] = getvalue(mp[:x][i,j]);
                end
                for par in 1:length(divSet[i])
                    yhat[i,par] = getvalue(mp[:y][i,par]);
                end
            end
            for ω in Ω
                θhat[ω] = getvalue(mp[:θ][ω]);
            end
            lbCost = getobjectivevalue(mp);
            push!(lbCostList,lbCost);
            push!(tcoreList,that);
            push!(xcoreList,xhat);
            push!(ycoreList,yhat);

            # examine how many cuts are tight at this solution, update the cutSel
            # cutSel,cutyn = examineCuts_count(disData,Ω,cutSel,cutSet,that,xhat,θhat,yhat,cutThreshold);
            # push!(cutynRec,length(cutyn));

            # update the master with new cutyn
            #yLim = limYselection(pData,H,that,divSet,radius);

            # generate cuts
            πdict = Dict();
            λdict = Dict();
            γdict = Dict();
            vk = Dict();
            θInt = Dict();
            ubTemp,θInt = ubCalP(pData,disData,Ω,xhat,that,Tmax1,1);
            if ubCost > ubTemp
                ubCost = ubTemp;
                tbest = copy(that);
                xbest = copy(xhat);
            end
            push!(ubCostList,ubCost);
            dataList = pmap(ω -> sub_divTDual(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict), Ω);

            # added the current relaxation
            mp1 = copy(mp);
            solve(mp1,relaxation = true);
            tcore = Dict();
            xcore = Dict();
            ycore = Dict();
            ycore1 = Dict();
            for i in pData.II
                tcore[i] = getvalue(mp1[:t][i]);
                for j in pData.Ji[i]
                    xcore[i,j] = getvalue(mp1[:x][i,j]);
                end
                for par in 1:length(divSet[i])
                    if (tcore[i] >= H[divSet[i][par].startH])&(tcore[i] < H[divSet[i][par].endH])
                        ycore1[i,par] = 1;
                    else
                        ycore1[i,par] = 0;
                    end
                    ycore[i,par] = getvalue(mp1[:y][i,par]);
                end
            end
            #dataList1 = pmap(ω -> sub_divTDualT(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore1), Ω);
            #dataList1 = pmap(ω -> sub_divTDualT2(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore,0.00001), Ω);
            dataList1 = pmap(ω -> sub_divTDualT3(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict,tcoreList,xcoreList,ycoreList), Ω);
            for ω in Ω
                πdict[ω] = dataList1[ω][1];
                λdict[ω] = dataList1[ω][2];
                γdict[ω] = dataList1[ω][3];
                vk[ω] = dataList1[ω][4];
            end
            ωTightCounter = 0;
            cutDual = [];
            for ω in Ω
                if vk[ω] - θhat[ω] > 1e-4*θhat[ω]
                    push!(cutDual,[ω,vk[ω],πdict[ω],λdict[ω],γdict[ω]]);
                    mp = addtxyCut(pData,ω,mp,πdict[ω],λdict[ω],γdict[ω],vk[ω],that,xhat,yhat,divSet);
                else
                    ωTightCounter += 1;
                end
            end
            if ωTightCounter == length(Ω)
                mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax,distanceDict,allSucc,0,yLim,0,cutyn);
                mp.solver = GurobiSolver();
                solve(mp);
                that = Dict();
                xhat = Dict();
                θhat = Dict();
                yhat = Dict();
                for i in pData.II
                    that[i] = getvalue(mp[:t][i]);
                    for j in pData.Ji[i]
                        xhat[i,j] = getvalue(mp[:x][i,j]);
                    end
                    for par in 1:length(divSet[i])
                        yhat[i,par] = getvalue(mp[:y][i,par]);
                    end
                end
                for ω in Ω
                    θhat[ω] = getvalue(mp[:θ][ω]);
                end
                lbCost = getobjectivevalue(mp);
                push!(lbCostList,lbCost);
                # generate cuts
                πdict = Dict();
                λdict = Dict();
                γdict = Dict();
                vk = Dict();
                θInt = Dict();
                ubTemp,θInt = ubCalP(pData,disData,Ω,xhat,that,Tmax1,1);
                if ubCost > ubTemp
                    ubCost = ubTemp;
                    tbest = copy(that);
                    xbest = copy(xhat);
                end
                push!(ubCostList,ubCost);
                dataList = pmap(ω -> sub_divT(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict), Ω);
                for ω in Ω
                    πdict[ω] = dataList[ω][1];
                    λdict[ω] = dataList[ω][2];
                    γdict[ω] = dataList[ω][3];
                    vk[ω] = dataList[ω][4];
                end
                ωTightCounter = 0;
                cutDual = [];
                for ω in Ω
                    if vk[ω] - θhat[ω] > 1e-4*θhat[ω]
                        push!(cutDual,[ω,vk[ω],πdict[ω],λdict[ω],γdict[ω]]);
                        mp = addtxyCut(pData,ω,mp,πdict[ω],λdict[ω],γdict[ω],vk[ω],that,xhat,yhat,divSet);
                    else
                        ωTightCounter += 1;
                    end
                end
                if ωTightCounter == length(Ω)
                    keepIter = false;
                    for i in pData.II
                        tlb[i] = that[i];
                        for j in pData.Ji[i]
                            xlb[i,j] = xhat[i,j];
                        end
                        for par in 1:length(divSet[i])
                            ylb[i,par] = yhat[i,par];
                        end
                    end
                    for ω in Ω
                        θlb[ω] = θhat[ω];
                    end
                else
                    cutSel,cutyn = examineCuts_count(disData,Ω,cutSel,cutSet,that,xhat,θhat,yhat,cutThreshold);
                    push!(cutynRec,length(cutyn));
                    # update the master with new cutyn
                    yLim = limYselection(pData,H,that,divSet,radius);
                    mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax,distanceDict,allSucc,1,yLim,0,cutyn);
                end
            else
                thatt = copy(that);
                xhatt = copy(xhat);
                yhatt = copy(yhat);
                push!(cutSet,[[thatt,xhatt,yhatt,divSet],cutDual]);
                # for l in 1:length(cutDual)
                #     cutSel[length(cutSet),l] = 0;
                # end
                mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax,distanceDict,allSucc,1,yLim,0,cutyn);
            end
        end
        tIter = toc();

        # need to come up with a rule to partition: gradient descent like binary search
        # check θInt vs. θhat: why the lower bound and the upper bound do not converge quickly --->
        # use the sub problem solution G to learn the b&b
        # also need to think up a way to tightening the cuts for each partition
        GFrac = Dict();
        for i in pData.II
            GFraciList = [ω for ω in Ω if (dataList[ω][5][i] < 1 - 1e-6)&(dataList[ω][5][i] > 1e-6)];
            if GFraciList != []
                GFrac[i] = [minimum(GFraciList),maximum(GFraciList)];
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

        # create new partition
        # if partOption == 1
        #     newPartition = [];
        #     for i in pData.II
        #         if GFrac[i] != []
        #             newItem = (i,Int(floor((GFrac[i][1] + GFrac[i][2])/2)));
        #             push!(newPartition,newItem);
        #         end
        #     end
        #     divSet,divDet = splitPar(divSet,divDet,newPartition);
        # elseif partOption == 2
        #     newPartition = [];
        #     for i in pData.II
        #         # obtain the CI by tRec
        #         meanTrec = mean(tRec[i]);
        #         CIu = meanTrec + 3*sqrt(sum((tRec[i][j] - meanTrec)^2 for j in 1:length(tRec[i]))/length(tRec[i]));
        #         CIuSet = [ω for ω in keys(H) if H[ω] > CIu];
        #         if CIuSet != []
        #             CIuH = minimum(CIuSet);
        #         else
        #             CIuH = maximum(keys(H));
        #         end
        #         CIl = meanTrec - 3*sqrt(sum((tRec[i][j] - meanTrec)^2 for j in 1:length(tRec[i]))/length(tRec[i]));
        #         CIlSet = [ω for ω in keys(H) if H[ω] < CIl];
        #         if CIlSet != []
        #             CIlH = maximum(CIlSet);
        #         else
        #             CIlH = 0;
        #         end
        #         newItem = (i,CIlH,CIuH);
        #         push!(newPartition,newItem);
        #     end
        #     divSet,divDet = splitPar_CI(divSet,divDet,newPartition);
        # end
        yLim = limYselection(pData,H,that,divSet,radius);
        cutThreshold += 5;
    end

    return tbest,xbest,lbCost,ubCost;
end
