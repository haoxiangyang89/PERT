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
    ubList,tHList,ubInc,tbest,xbest,θbest,textList,xextList = iniPart(pData,disData,Ω,sN,MM,1,noThreads);
    lbCost = -Inf;
    lbCostList = [];
    global ubCost = ubInc;
    global ubCostList = [ubCost];

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

    lbHist = [];
    ubHist = [];
    timeHist = [];
    cutHist = [];
    intSolHist = [];
    yhistList = [];

    # set up a tree list
    treeList = [];
    push!(treeList,[divSet,divDet,lbCost,cutSet,-1]);
    global batchNo = 5;
    global lbOverAll = -Inf;


    function runPara(pData,disData,treeList,tcoreList,xcoreList,ycoreList,H,Ω,lDict,allSucc,distanceDict,ubCost,tbest,xbest,batchNo)
        npList = [ib for ib in 2:batchNo+1];
        global keepIter = true;
        global noTh = div(noThreads,batchNo);
        @sync begin
            for p in npList
                @async begin
                    while true
                        # if all nodes are processed and no nodes are being processed, exit
                        boolFinished = true;
                        for l in 1:length(treeList)
                            if treeList[l][5] != 1
                                boolFinished = false;
                            end
                        end
                        println(boolFinished," ",keepIter," ",[treeList[l][5] for l in 1:length(treeList)]);
                        if (boolFinished) || (!(keepIter))
                            println("break");
                            break
                        end
                        openNodes = [(treeList[l][3],l) for l in 1:length(treeList) if treeList[l][5] == -1];
                        if openNodes != []
                            selectNode = sort(openNodes, by = x -> x[1])[1][2];
                            println(p," ",selectNode);
                            treeList[selectNode][5] = 0;
                            mpSolveInfo = remotecall_fetch(solveMP_para,p,[pData,disData,H,Ω,lDict,allSucc,distanceDict,treeList[selectNode][1],
                                                        treeList[selectNode][2],treeList[selectNode][3],tcoreList,xcoreList,ycoreList,
                                                        ubCost,tbest,xbest,noTh]);
                            println(mpSolveInfo);
                            treeList[selectNode][5] = 1;
                            for ib in 1:length(divSetList)
                                if mpSolveInfo[ib][6] < ubCost
                                    ubCost = mpSolveInfo[ib][6];
                                    tbest = mpSolveInfo[ib][4];
                                    xbest = mpSolveInfo[ib][5];
                                end
                                append!(tcoreList,mpSolveInfo[ib][8]);
                                append!(xcoreList,mpSolveInfo[ib][9]);
                                append!(ycoreList,mpSolveInfo[ib][10]);
                            end
                            # compare the current lb with the current best lb
                            lbOverAll = minimum([treeList[l][3] for l in 1:length(treeList) if treeList[l][5] != 1]);
                            mpStatus = mpSolveInfo[ibatch][1];
                            mpObj = mpSolveInfo[ibatch][2];
                            if mpStatus == :Optimal
                                if mpObj < lbOverAll
                                    lbOverAll = mpObj;
                                end
                            end
                            if (ubCost - lbOverAll)/ubCost < ϵ
                                keepIter = false;
                            else
                                if (mpSolveInfo[ibatch][2] < ubCost)&(mpSolveInfo[ibatch][1] == :Optimal)
                                    # branch the current node
                                    GList = mpSolveInfo[ibatch][3];
                                    cutSet = mpSolveInfo[ibatch][7];
                                    GCurrent = GList[length(GList)];
                                    GFrac = Dict();
                                    lbCost = mpSolveInfo[ibatch][2];
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
                                    divSet1,divDet1,divSet2,divDet2 = breakDiv(pData,disData,H,divSetList[ibatch],divDetList[ibatch],lGFracInd,locBreak,distanceDict);
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
                                    push!(treeList,[divSet1,divDet1,lbCost,deepcopy(cutSet),0]);

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
                                    push!(treeList,[divSet2,divDet2,lbCost,deepcopy(cutSet),0]);
                                end
                            end
                        else
                            # sleep if there is no current available nodes
                            println(p);
                            remotecall_fetch(sleep, p, 10);
                        end
                    end
                end
            end
        end
    end
    #     while length([iTree for iTree in 1:length(treeList) if treeList[iTree][5] == 0]) > 0
    #         for bn in availableWorker
    #             @async begin
    #                     if length([iTree for iTree in 1:length(treeList) if treeList[iTree][5] == 0]) > 0
    #                         minTree = Inf;
    #                         minTreeInd = -1;
    #                         for iTree in 1:length(treeList)
    #                             if (treeList[iTree][3] < minTree)&(treeList[iTree][5] == 0)
    #                                 minTree = treeList[iTree][3];
    #                                 minTreeInd = iTree;
    #                             end
    #                         end
    #                         if minTreeInd != -1
    #                             push!(divSetList,deepcopy(treeList[minTreeInd][1]));
    #                             push!(divDetList,deepcopy(treeList[minTreeInd][2]));
    #                             push!(cutSetList,deepcopy(treeList[minTreeInd][4]));
    #
    #                             treeList[minTreeInd][5] = 1;
    #                         end
    #                         if length([iTree for iTree in 1:length(treeList) if treeList[iTree][5] == 0]) == 0
    #                             lbOverAll = Inf;
    #                         else
    #                             lbOverAll = minimum([treeList[iTree][3] for iTree in 1:length(treeList) if treeList[iTree][5] == 0]);
    #                         end
    #                         mpSolveList = solveMP_para(pData,disData,H,Ω,lDict,allSucc,distanceDict,
    #                             divSetList[ib],divDetList[ib],cutSetList[ib],tcoreList,xcoreList,ycoreList,
    #                             ubCost,tbest,xbest,noTh);
    #                         # record the results from solving the node
    #                         for ib in 1:length(divSetList)
    #                             if mpSolveList[ib][6] < ubCost
    #                                 ubCost = mpSolveList[ib][6];
    #                                 tbest = mpSolveList[ib][4];
    #                                 xbest = mpSolveList[ib][5];
    #                             end
    #                             append!(tcoreList,mpSolveList[ib][8]);
    #                             append!(xcoreList,mpSolveList[ib][9]);
    #                             append!(ycoreList,mpSolveList[ib][10]);
    #                         end
    #                         push!(ubHist,ubCost);
    #                         for ibatch in 1:length(divSetList)
    #                             mpStatus = mpSolveList[ibatch][1];
    #                             mpObj = mpSolveList[ibatch][2];
    #                             if mpStatus == :Optimal
    #                                 if mpObj < lbOverAll
    #                                     lbOverAll = mpObj;
    #                                 end
    #                             end
    #                         end
    #                         if (ubCost - lbOverAll)/ubCost < ϵ
    #                             keepIter = false;
    #                         else
    #                         #        cutSel = examineCuts_count_2(disData,Ω,cutSet,divSet,tCurrent,xCurrent,θCurrent,yCurrent);
    #                         #        cutSetNew = selectCuts2(cutSet,cutSel);
    #                                 # if the current lower bound is larger than the upper bound, prune the node
    #                             for ibatch in 1:length(divSetList)
    #                                 if (mpSolveList[ibatch][2] < ubCost)&(mpSolveList[ibatch][1] == :Optimal)
    #                                     GList = mpSolveList[ibatch][3];
    #                                     cutSet = mpSolveList[ibatch][7];
    #                                     GCurrent = GList[length(GList)];
    #                                     GFrac = Dict();
    #                                     lbCost = mpSolveList[ibatch][2];
    #                                     for i in pData.II
    #                                         GFraciList = [disData[ω].H for ω in Ω if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)];
    #                                         if GFraciList != []
    #                                             GFrac[i] = [HRev[minimum(GFraciList)],HRev[maximum(GFraciList)]];
    #                                         else
    #                                             GFrac[i] = [];
    #                                         end
    #                                     end
    #                                     newPartition = [];
    #                                     for i in pData.II
    #                                         if GFrac[i] != []
    #                                             newItem = (i,GFrac[i][1],GFrac[i][2]);
    #                                             push!(newPartition,newItem);
    #                                         end
    #                                     end
    #                                     # select the largest GFrac coverage
    #                                     largestGFrac = -Inf;
    #                                     lGFracInd = -1;
    #                                     for i in pData.II
    #                                         if GFrac[i] != []
    #                                             if GFrac[i][2] - GFrac[i][1] > largestGFrac
    #                                                 largestGFrac = GFrac[i][2] - GFrac[i][1];
    #                                                 lGFracInd = i;
    #                                             end
    #                                         end
    #                                     end
    #                                     locBreak = Int64(floor((GFrac[lGFracInd][1]*2/3 + GFrac[lGFracInd][2]*1/3)));
    #                                     divSet1,divDet1,divSet2,divDet2 = breakDiv(pData,disData,H,divSetList[ibatch],divDetList[ibatch],lGFracInd,locBreak,distanceDict);
    #                                     lGFracInd1 = -1;
    #                                     largest1 = -Inf;
    #                                     fracTopLG = -1;
    #                                     fracBotLG = -1;
    #                                     for i in pData.II
    #                                         if (!(i in allSucc[lGFracInd]))&(!(lGFracInd in allSucc[i]))&(i != lGFracInd)
    #                                             if GFrac[i] != []
    #                                                 fracTop = min(GFrac[i][2],maximum([divSet1[i][par].endH for par in 1:length(divSet1[i]) if divDet1[i][par] == 0]));
    #                                                 fracBot = max(GFrac[i][1],minimum([divSet1[i][par].startH for par in 1:length(divSet1[i]) if divDet1[i][par] == 0]));
    #                                                 #println(i," ",fracTop," ",fracBot);
    #                                                 if fracTop - fracBot > largest1
    #                                                     largest1 = fracTop - fracBot;
    #                                                     lGFracInd1 = i;
    #                                                     fracTopLG = fracTop;
    #                                                     fracBotLG = fracBot;
    #                                                 end
    #                                             end
    #                                         end
    #                                     end
    #                                     newPartition1 = [(lGFracInd1,fracBotLG,fracTopLG)];
    #                                     divSet1,divDet1 = splitPar3(divSet1,divDet1,newPartition1);
    #                                     push!(treeList,[divSet1,divDet1,lbCost,deepcopy(cutSet),0]);
    #
    #                                     lGFracInd2 = -1;
    #                                     largest2 = -Inf;
    #                                     fracTopLG = -1;
    #                                     fracBotLG = -1;
    #                                     for i in pData.II
    #                                         if (!(i in allSucc[lGFracInd]))&(!(lGFracInd in allSucc[i]))&(i != lGFracInd)
    #                                             if GFrac[i] != []
    #                                                 fracTop = min(GFrac[i][2],maximum([divSet2[i][par].endH for par in 1:length(divSet2[i]) if divDet2[i][par] == 0]));
    #                                                 fracBot = max(GFrac[i][1],minimum([divSet2[i][par].startH for par in 1:length(divSet2[i]) if divDet2[i][par] == 0]));
    #                                                 if fracTop - fracBot > largest2
    #                                                     largest2 = fracTop - fracBot;
    #                                                     lGFracInd2 = i;
    #                                                     fracTopLG = fracTop;
    #                                                     fracBotLG = fracBot;
    #                                                 end
    #                                             end
    #                                         end
    #                                     end
    #                                     newPartition2 = [(lGFracInd2,fracBotLG,fracTopLG)];
    #                                     divSet2,divDet2 = splitPar3(divSet2,divDet2,newPartition2);
    #                                     push!(treeList,[divSet2,divDet2,lbCost,deepcopy(cutSet),0]);
    #                                 end
    #                             end
    #                         end
    #                     end
    #                 end
    #         end
    #     end
    # end

    while keepIter
        # select the node with lowest lower bound
        divSetList = [];
        divDetList = [];
        cutSetList = [];
        for ibatch in 1:batchNo
            minTree = Inf;
            minTreeInd = -1;
            for iTree in 1:length(treeList)
                if (treeList[iTree][3] < minTree)&(treeList[iTree][5] == 0)
                    minTree = treeList[iTree][3];
                    minTreeInd = iTree;
                end
            end
            if minTreeInd != -1
                push!(divSetList,deepcopy(treeList[minTreeInd][1]));
                push!(divDetList,deepcopy(treeList[minTreeInd][2]));
                push!(cutSetList,deepcopy(treeList[minTreeInd][4]));

                treeList[minTreeInd][5] = 1;
            end
        end
        # if there are open nodes, record the smallest cost, else record +inf
        if [treeList[l][3] for l in 1:length(treeList) if treeList[l][5] == 0] != []
            lbOverAll = minimum([treeList[l][3] for l in 1:length(treeList) if treeList[l][5] == 0]);
        else
            lbOverAll = Inf;
        end

        noTh = div(noThreads,length(divSetList));
        tic();
        #mpStatus = solve(mp);
        # mpSolveList = asyncmap(ib -> solveMP_para(pData,disData,H,Ω,lDict,allSucc,distanceDict,
        #     divSetList[ib],divDetList[ib],cutSetList[ib],tcoreList,xcoreList,ycoreList,
        #     ubCost,tbest,xbest,noTh),1:length(divSetList); ntasks = batchNo);
        mpSolveList = pmap(ib -> solveMP_para(pData,disData,H,Ω,lDict,allSucc,distanceDict,
            divSetList[ib],divDetList[ib],cutSetList[ib],tcoreList,xcoreList,ycoreList,
            ubCost,tbest,xbest,noTh),1:length(divSetList)); # some of the function is not run with full number of cores
        tIter = toc();
        push!(timeHist,tIter);
        for ib in 1:length(divSetList)
            if mpSolveList[ib][6] < ubCost
                ubCost = mpSolveList[ib][6];
                tbest = mpSolveList[ib][4];
                xbest = mpSolveList[ib][5];
            end
            append!(tcoreList,mpSolveList[ib][8]);
            append!(xcoreList,mpSolveList[ib][9]);
            append!(ycoreList,mpSolveList[ib][10]);
        end
        push!(ubHist,ubCost);
        for ibatch in 1:length(divSetList)
            mpStatus = mpSolveList[ibatch][1];
            mpObj = mpSolveList[ibatch][2];
            if mpStatus == :Optimal
                if mpObj < lbOverAll
                    lbOverAll = mpObj;
                end
            end
        end

            # need to come up with a rule to partition: gradient descent like binary search
            # check θInt vs. θhat: why the lower bound and the upper bound do not converge quickly --->
            # use the sub problem solution G to learn the b&b
            # also need to think up a way to tightening the cuts for each partition
            #push!(intSolHist,length(ubCostList));
        if (ubCost - lbOverAll)/ubCost < ϵ
            keepIter = false;
        else
        #        cutSel = examineCuts_count_2(disData,Ω,cutSet,divSet,tCurrent,xCurrent,θCurrent,yCurrent);
        #        cutSetNew = selectCuts2(cutSet,cutSel);
                # if the current lower bound is larger than the upper bound, prune the node
            for ibatch in 1:length(divSetList)
                if (mpSolveList[ibatch][2] < ubCost)&(mpSolveList[ibatch][1] == :Optimal)
                    GList = mpSolveList[ibatch][3];
                    cutSet = mpSolveList[ibatch][7];
                    GCurrent = GList[length(GList)];
                    GFrac = Dict();
                    lbCost = mpSolveList[ibatch][2];
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
                    divSet1,divDet1,divSet2,divDet2 = breakDiv(pData,disData,H,divSetList[ibatch],divDetList[ibatch],lGFracInd,locBreak,distanceDict);
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
                    #newPartition1 = [(lGFracInd,locBreak,GFrac[lGFracInd][2])];
                    divSet1,divDet1 = splitPar3(divSet1,divDet1,newPartition1);
                    push!(treeList,[divSet1,divDet1,lbCost,deepcopy(cutSet),0]);

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
                    #newPartition2 = [(lGFracInd,GFrac[lGFracInd][1],locBreak)];
                    divSet2,divDet2 = splitPar3(divSet2,divDet2,newPartition2);
                    push!(treeList,[divSet2,divDet2,lbCost,deepcopy(cutSet),0]);
                end
                # divSet,divDet = splitPar(divSet,divDet,newPartition);
        #        push!(cutHist,sum(length(cutSet[l][2]) for l in 1:length(cutSet)));
        #        cutSet = deepcopy(cutSetNew);
            end
        end
    end

    # need a cut selection process within the callback
    return tbest,xbest,ubCost,lbOverAll;
end
