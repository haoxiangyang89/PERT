# partition by t
function revPar(prevPartSet,currentPart)
    # i specific
    corrPart = 0;
    for par in 1:length(prevPartSet)
        if (currentPart.startH >= prevPartSet[par].startH)&(currentPart.endH <= prevPartSet[par].endH)
            corrPart = par;
        end
    end
    return corrPart;
end

function splitPar(PartSet,PartDet,splitInfo)
    # i generic
    newPartSet = copy(PartSet);
    newPartDet = copy(PartDet);
    for (i,splitStart,splitEnd) in splitInfo
        partSetiTemp = [];
        partDetiTemp = [];
        splitPos = Int(floor((splitStart + splitEnd)/2));
        for par in 1:length(PartSet[i])
            # split the current partition
            if (splitEnd <= PartSet[i][par].endH)&(splitStart >= PartSet[i][par].startH)
                part1 = partType(PartSet[i][par].startH,splitPos);
                part2 = partType(splitPos,PartSet[i][par].endH);
                push!(partSetiTemp,part1);
                push!(partSetiTemp,part2);
                push!(partDetiTemp,0);
                push!(partDetiTemp,0);
            else
                push!(partSetiTemp,PartSet[i][par]);
                push!(partDetiTemp,PartDet[i][par]);
            end
        end
        newPartSet[i] = partSetiTemp;
        newPartDet[i] = partDetiTemp;
    end
    return newPartSet,newPartDet;
end

function splitPar3(PartSet,PartDet,splitInfo)
    # i generic
    newPartSet = copy(PartSet);
    newPartDet = copy(PartDet);
    for (i,splitStart,splitEnd) in splitInfo
        partSetiTemp = [];
        partDetiTemp = [];
        for par in 1:length(PartSet[i])
            # split the current partition
            if (splitEnd <= PartSet[i][par].endH)&(splitStart >= PartSet[i][par].startH)
                mod3 = mod(PartSet[i][par].endH - PartSet[i][par].startH,3);
                len3 = div(PartSet[i][par].endH - PartSet[i][par].startH,3);
                if len3 >= 1
                    # if the partition is larger than 3 items in between
                    if mod3 == 0
                        part1 = partType(PartSet[i][par].startH,PartSet[i][par].startH + len3);
                        part2 = partType(PartSet[i][par].startH + len3,PartSet[i][par].startH + len3*2);
                        part3 = partType(PartSet[i][par].startH + len3*2,PartSet[i][par].endH);
                    elseif mod3 == 1
                        part1 = partType(PartSet[i][par].startH,PartSet[i][par].startH + len3);
                        part2 = partType(PartSet[i][par].startH + len3,PartSet[i][par].startH + len3*2);
                        part3 = partType(PartSet[i][par].startH + len3*2,PartSet[i][par].endH);
                    else
                        part1 = partType(PartSet[i][par].startH,PartSet[i][par].startH + len3);
                        part2 = partType(PartSet[i][par].startH + len3,PartSet[i][par].startH + len3*2 + 1);
                        part3 = partType(PartSet[i][par].startH + len3*2 + 1,PartSet[i][par].endH);
                    end
                    push!(partSetiTemp,part1);
                    push!(partSetiTemp,part2);
                    push!(partSetiTemp,part3);
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
                else
                    splitPos = Int(floor((splitStart + splitEnd)/2));
                    part1 = partType(PartSet[i][par].startH,splitPos);
                    part2 = partType(splitPos,PartSet[i][par].endH);
                    push!(partSetiTemp,part1);
                    push!(partSetiTemp,part2);
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
                end
            else
                push!(partSetiTemp,PartSet[i][par]);
                push!(partDetiTemp,PartDet[i][par]);
            end
        end
        newPartSet[i] = partSetiTemp;
        newPartDet[i] = partDetiTemp;
    end
    return newPartSet,newPartDet;
end

function revisePar(pData,disData,PartSet,PartDet,ubInfo,lbInfo)
    H = Dict();
    H[0] = 0;
    H[length(Ω)+1] = Tmax;
    for ω in Ω
        H[ω] = disData[ω].H;
    end

    newPartSet = copy(PartSet);
    newPartDet = copy(PartDet);
    # for each activity
    for i in pData.II
        partSetiTemp = [];
        partDetiTemp = [];
        for par in 1:length(PartSet[i])
            if PartDet[i][par] == 0
                set1 = [ω for ω in (PartSet[i][par].startH + 1):(PartSet[i][par].endH - 1) if H[ω] < lbInfo[i]];
                set3 = [ω for ω in (PartSet[i][par].startH + 1):(PartSet[i][par].endH - 1) if H[ω] >= ubInfo[i]];
                if set1 != []
                    # split the current set into parts with the first det as 1
                    set1start = PartSet[i][par].startH;
                    set1end = maximum(set1);
                    if set1end == PartSet[i][par].endH
                        push!(partSetiTemp,PartSet[i][par]);
                        push!(partDetiTemp,1);
                    else
                        push!(partSetiTemp,partType(set1start,set1end));
                        push!(partDetiTemp,1);
                        set2start = set1end;
                        if set3 == []
                            set2end = PartSet[i][par].endH;
                            push!(partSetiTemp,partType(set2start,set2end));
                            push!(partDetiTemp,0);
                        else
                            set2end = minimum(set3);
                            push!(partSetiTemp,partType(set2start,set2end));
                            push!(partDetiTemp,0);
                            set3start = set2end;
                            set3end = PartSet[i][par].endH;
                            push!(partSetiTemp,partType(set3start,set3end));
                            push!(partDetiTemp,-1);
                        end
                    end
                    # split the current set into parts with the first det as 0
                else
                    set1start = PartSet[i][par].startH;
                    if set3 == []
                        push!(partSetiTemp,PartSet[i][par]);
                        push!(partDetiTemp,0);
                    else
                        set1end = minimum(set3);
                        push!(partSetiTemp,partType(set1start,set1end));
                        push!(partDetiTemp,0);
                        set3start = set1end;
                        set3end = PartSet[i][par].endH;
                        push!(partSetiTemp,partType(set3start,set3end));
                        push!(partDetiTemp,-1);
                    end
                end
            else
                push!(partSetiTemp,PartSet[i][par]);
                push!(partDetiTemp,PartDet[i][par]);
            end
        end
        newPartSet[i] = partSetiTemp;
        newPartDet[i] = partDetiTemp;
    end
    return newPartSet,newPartDet;
end

function splitPar_CI(PartSet,PartDet,splitInfo,forceZero = false)
    # i generic
    newPartSet = copy(PartSet);
    newPartDet = copy(PartDet);
    for (i,splitStart,splitEnd) in splitInfo
        partSetiTemp = [];
        partDetiTemp = [];
        for par in 1:length(PartSet[i])
            # split the current partition, essentially split twice using the start and the end
            if (splitStart < PartSet[i][par].endH)&(splitStart > PartSet[i][par].startH)&(PartDet[i][par] == 0)
                if (splitEnd < PartSet[i][par].endH)&(splitEnd > PartSet[i][par].startH)
                    part1 = partType(PartSet[i][par].startH,splitStart);
                    part2 = partType(splitStart,splitEnd);
                    part3 = partType(splitEnd,PartSet[i][par].endH);
                    push!(partSetiTemp,part1);
                    push!(partSetiTemp,part2);
                    push!(partSetiTemp,part3);
                    if forceZero
                        push!(partDetiTemp,1);
                        push!(partDetiTemp,0);
                        push!(partDetiTemp,-1);
                    else
                        push!(partDetiTemp,0);
                        push!(partDetiTemp,0);
                        push!(partDetiTemp,0);
                    end
                else
                    part1 = partType(PartSet[i][par].startH,splitStart);
                    part2 = partType(splitStart,PartSet[i][par].endH);
                    push!(partSetiTemp,part1);
                    push!(partSetiTemp,part2);
                    if forceZero
                        push!(partDetiTemp,1);
                        push!(partDetiTemp,0);
                    else
                        push!(partDetiTemp,0);
                        push!(partDetiTemp,0);
                    end
                end
            else
                if (splitEnd < PartSet[i][par].endH)&(splitEnd > PartSet[i][par].startH)
                    part1 = partType(PartSet[i][par].startH,splitEnd);
                    part2 = partType(splitEnd,PartSet[i][par].endH);
                    push!(partSetiTemp,part1);
                    push!(partSetiTemp,part2);
                    if forceZero
                        push!(partDetiTemp,0);
                        push!(partDetiTemp,-1);
                    else
                        push!(partDetiTemp,0);
                        push!(partDetiTemp,0);
                    end
                else
                    push!(partSetiTemp,PartSet[i][par]);
                    push!(partDetiTemp,PartDet[i][par]);
                end
            end
        end
        newPartSet[i] = partSetiTemp;
        newPartDet[i] = partDetiTemp;
    end
    return newPartSet,newPartDet;
end

function prunePart(pData,disData,Ω,divSet,divDet,cutSet,Tmax,distanceDict,allSucc,currentUB)
    # try to solve the linear relaxation of the master to prune some partition
    for i in pData.II
        for par in 1:length(divSet[i])
            if divDet[i][par] == 0
                mpTemp = createMaster_DivRel(pData,disData,Ω,divSet,divDet,cutSet,Tmax,distanceDict,allSucc,[i,par]);
                #mpTemp.solver = GurobiSolver();
                mpStatus = solve(mpTemp);
                println(i," ",par," ",getobjectivevalue(mpTemp));
                if mpStatus == :Optimal
                    if getobjectivevalue(mpTemp) > currentUB
                        divDet[i][par] = -1;
                    end
                else
                    divDet[i][par] = -1;
                end
            end
        end
    end
    return divDet;
end

# find all subsequent nodes of activity i
function findSuccAll(pData)
    allSucc = Dict();
    for i in pData.II
        succList = copy(pData.Succ[i]);
        succAll = [];
        while succList != []
            currentSucc = succList[1];
            push!(succAll,currentSucc);
            shift!(succList);
            for j in pData.Succ[currentSucc]
                if !((j in succAll) | (j in succList))
                    push!(succList,j);
                end
            end
        end
        allSucc[i] = succAll;
    end
    return allSucc;
end

function iniPart(pData,disData,Ω,sN,MM,returnOpt = 0,noThreads = 30)
    Tmax = disData[length(Ω)].H + longestPath(pData)[0];
    H = Dict();
    H[0] = 0;
    H[length(Ω)+1] = Tmax;
    for ω in Ω
        H[ω] = disData[ω].H;
    end

    # sample sN scenarios and solve the small extensive formulation
    tList = Dict();
    for i in pData.II
        tList[i] = [];
    end
    ubList = [];
    ubMin = Inf;
    tBest = Dict();
    xBest = Dict();
    θBest = Dict();
    textList = [];
    xextList = [];
    m = 1;
    while m <= MM
        try
            #randList = rand(Ω,sN);
            disData1 = Dict();
            Ω1 = 1:sN;
            randList = [(ω - 1)*MM + m for ω in Ω1];
            for i in Ω1
                disData1[i] = deepcopy(disData[randList[i]]);
                disData1[i].prDis = (1 - pData.p0)/length(Ω1);
            end
            text,xext,fext,gext,mext = extForm_cheat(pData,disData1,Ω1,1e-4,999999,noThreads);
            for i in pData.II
                if abs(text[i]) <= 1e-5
                    text[i] = 0.0;
                end
            end
            push!(textList,text);
            push!(xextList,xext);
            ubext,cωList = ubCalP(pData,disData,Ω,xext,text,999999,1);
            push!(ubList,ubext);
            # record the best solution
            if ubext < ubMin
                ubMin = ubext;
                for i in pData.II
                    tBest[i] = text[i];
                    for j in pData.Ji[i]
                        xBest[i,j] = xext[i,j];
                    end
                end
                for ω in Ω
                    θBest[ω] = cωList[ω];
                end
            end
            for i in pData.II
                push!(tList[i],text[i]);
            end
            m += 1;
        catch
            sleep(30);
        end
    end

    # return the min-max decomposition
    tHList = [];
    for i in pData.II
        tHstart = maximum([ω for ω in keys(H) if H[ω] <= minimum(tList[i])]);
        tHend = minimum([ω for ω in keys(H) if H[ω] >= maximum(tList[i])]);
        push!(tHList,[i,tHstart,tHend]);
    end
    if returnOpt == 0
        return ubList,tHList,ubMin,tBest,xBest,θBest;
    else
        return ubList,tHList,ubMin,tBest,xBest,θBest,textList,xextList;
    end
end

function splitAny(PartSet,PartDet,Ω,breakPoints)
    # breakPoints is a dictionary with keys from pData.II
    newPartSet = deepcopy(PartSet);
    newPartDet = deepcopy(PartDet);
    for i in keys(breakPoints)
        partSetiTemp = [];
        partDetiTemp = [];
        for par in 1:length(PartSet[i])
            # split the current partition, essentially split twice using the start and the end
            bps = sort([j for j in breakPoints[i] if (j < PartSet[i][par].endH)&(j > PartSet[i][par].startH)]);
            if bps != []
                currentStart = PartSet[i][par].startH;
                for j in bps
                    partCurrent = partType(currentStart,j);
                    push!(partSetiTemp,partCurrent);
                    push!(partDetiTemp,0);
                    currentStart = j;
                end
                currentEnd = PartSet[i][par].endH;
                partCurrent = partType(currentStart,currentEnd);
                push!(partSetiTemp,partCurrent);
                push!(partDetiTemp,0);
            else
                push!(partSetiTemp,PartSet[i][par]);
                push!(partDetiTemp,PartDet[i][par]);
            end
        end
        newPartSet[i] = partSetiTemp;
        newPartDet[i] = partDetiTemp;
    end
    return newPartSet,newPartDet;
end

function splitPrep3(pData,disData,Ω,H,HRev,GCurrent,divSet,divDet)
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
    divSetNew,divDetNew = splitPar3(divSet,divDet,newPartition);
    return divSetNew,divDetNew;
end

function splitPrepld(pData,disData,Ω,H,GCurrent,divSet,divDet,θlp,θInt,nSplit)
    # split the partition according to the largest deviation between the lp relaxation and the ip sub
    θDiff = [θInt[ω] - θlp[ω] for ω in Ω];
    θDiffPerm = sortperm(θDiff,rev = true);
    # initialize the breakPoints dictionary
    divDict = Dict();
    for i in pData.II
        divDict[i] = [];
    end
    for n in 1:nSplit
        ωn = Ω[θDiffPerm[n]];
        for i in pData.II
            # G[i] is fractional
            if (GCurrent[ωn][i] < 1 - 1e-6)&(GCurrent[ωn][i] > 1e-6)
                push!(divDict[i],ωn);
            end
        end
    end
    divSetNew,divDetNew = splitAny(divSet,divDet,Ω,divDict);
    return divSetNew,divDetNew;
end

function splitPrepld2(pData,disData,Ω,H,GCurrent,tCurrent,divSet,divDet,θlp,θInt,nSplit)
    # split the partition according to the largest fraction impact (= G deviation * disData.d)
    divDict = Dict();
    for i in pData.II
        divDict[i] = [];
    end
    for i in pData.II
        θDiv = [];
        for par in 1:length(divSet[i])
            if divDet[i][par] == 0
                for ω in (divSet[i][par].startH + 1):(divSet[i][par].endH - 1)
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
    divSetNew,divDetNew = splitAny(divSet,divDet,Ω,divDict);
    return divSetNew,divDetNew;
end

function splitPrepSmart(pData,disData,Ω,H,GCurrent,tCurrent,divSet,divDet,θlp,θInt,nSplit)
    # select the partition according to both sparsity and scenario improvement
    divSetNew = deepcopy(divSet);
    divDetNew = deepcopy(divDet);
    divΩ = [];
    n = 1;
    contBool = true;
    while (n <= nSplit)&(contBool)
        divDict = Dict();
        for i in pData.II
            divDict[i] = [];
        end
        # select the split in a greedy manner
        θDiv = [];
        for ω in Ω
            if (θInt[ω] - θlp[ω] > θInt[ω]*1e-3)&(!(ω in divΩ))
                # if we select ω to branch on
                θgain = θInt[ω] - θlp[ω];
                θgainB = [];
                for i in pData.II
                    if i != 0
                        # obtain the original partition diameter
                        currentpar = -1;
                        for par in 1:length(divSetNew[i])
                            if (ω > divSetNew[i][par].startH) & (ω < divSetNew[i][par].endH)
                                currentpar = par;
                            end
                        end
                        if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)&(currentpar != -1)
                            if divDetNew[i][currentpar] == 0
                                # if G is fractional
                                oriD = H[divSetNew[i][currentpar].endH] - H[divSetNew[i][currentpar].startH];
                                for ωn in (divSetNew[i][currentpar].startH+1):(ω-1)
                                    currentD = H[ω] - H[divSetNew[i][currentpar].startH];
                                    if tCurrent[i] > H[ωn] + 1e-6
                                        push!(θgainB,(1 - GCurrent[ω][i])*disData[ωn].d[i]*(1/currentD - 1/oriD));
                                    elseif tCurrent[i] < H[ω] - 1e-6
                                        push!(θgainB,GCurrent[ω][i]*disData[ωn].d[i]*(1/currentD - 1/oriD));
                                    else
                                        push!(θgainB,min(1 - GCurrent[ω][i],GCurrent[ωn][i])*disData[ωn].d[i]*(1/currentD - 1/oriD));
                                    end
                                end
                                for ωn in (ω+1):(divSetNew[i][currentpar].endH-1)
                                    currentD = H[divSetNew[i][currentpar].endH] - H[ω];
                                    if tCurrent[i] > H[ωn] + 1e-6
                                        push!(θgainB,(1 - GCurrent[ωn][i])*disData[ωn].d[i]*(1/currentD - 1/oriD));
                                    elseif tCurrent[i] < H[ω] - 1e-6
                                        push!(θgainB,GCurrent[ωn][i]*disData[ωn].d[i]*(1/currentD - 1/oriD));
                                    else
                                        push!(θgainB,min(1 - GCurrent[ωn][i],GCurrent[ωn][i])*disData[ωn].d[i]*(1/currentD - 1/oriD));
                                    end
                                end
                            end
                        end
                    end
                end
                if θgainB != []
                    θgain += maximum(θgainB);
                end
                push!(θDiv,(θgain,ω));
            end
        end
        if θDiv != []
            estθDiff,ωSelect = maximum(θDiv);
            for i in pData.II
                # G[i] is fractional
                if (GCurrent[ωSelect][i] < 1 - 1e-6)&(GCurrent[ωSelect][i] > 1e-6)
                    push!(divDict[i],ωSelect);
                end
            end
            divSetNew,divDetNew = splitAny(divSetNew,divDetNew,Ω,divDict);
            push!(divΩ,ωSelect);
        else
            contBool = false;
        end
        n += 1;
    end
    return divSetNew,divDetNew;
end

function splitPrepSmart2(pData,disData,Ω,H,GCurrent,tCurrent,divSet,divDet,θlp,θInt,nSplit)
    # select the partition according to both sparsity and scenario improvement
    divSetNew = deepcopy(divSet);
    divDetNew = deepcopy(divDet);
    divΩ = [];
    n = 1;
    while (n <= nSplit)
        divDict = Dict();
        for i in pData.II
            divDict[i] = [];
        end
        # select the split in a greedy manner
        for i in pData.II
            θDiv = [];
            if i != 0
                for ω in Ω
                    if (θInt[ω] - θlp[ω] > θInt[ω]*1e-3)
                        # if we select ω to branch on
                        θgain = 0;
                        # obtain the original partition diameter
                        currentpar = -1;
                        for par in 1:length(divSetNew[i])
                            if (ω > divSetNew[i][par].startH) & (ω < divSetNew[i][par].endH)
                                currentpar = par;
                            end
                        end
                        if (GCurrent[ω][i] < 1 - 1e-6)&(GCurrent[ω][i] > 1e-6)&(currentpar != -1)
                            if divDetNew[i][currentpar] == 0
                                # if G is fractional
                                oriD = H[divSetNew[i][currentpar].endH] - H[divSetNew[i][currentpar].startH];
                                for ωn in (divSetNew[i][currentpar].startH+1):(ω-1)
                                    currentD = H[ω] - H[divSetNew[i][currentpar].startH];
                                    if tCurrent[i] > H[ωn] + 1e-6
                                        θgain += (1 - GCurrent[ω][i])*disData[ωn].d[i]*(1/currentD - 1/oriD);
                                    elseif tCurrent[i] < H[ω] - 1e-6
                                        θgain += GCurrent[ω][i]*disData[ωn].d[i]*(1/currentD - 1/oriD);
                                    else
                                        θgain += min(1 - GCurrent[ω][i],GCurrent[ωn][i])*disData[ωn].d[i]*(1/currentD - 1/oriD);
                                    end
                                end
                                for ωn in (ω+1):(divSetNew[i][currentpar].endH-1)
                                    currentD = H[divSetNew[i][currentpar].endH] - H[ω];
                                    if tCurrent[i] > H[ωn] + 1e-6
                                        θgain += (1 - GCurrent[ωn][i])*disData[ωn].d[i]*(1/currentD - 1/oriD);
                                    elseif tCurrent[i] < H[ω] - 1e-6
                                        θgain += GCurrent[ωn][i]*disData[ωn].d[i]*(1/currentD - 1/oriD);
                                    else
                                        θgain += min(1 - GCurrent[ωn][i],GCurrent[ωn][i])*disData[ωn].d[i]*(1/currentD - 1/oriD);
                                    end
                                end
                                push!(θDiv,(θgain,ω));
                            end
                        end
                    end
                end
                if θDiv != []
                    estθDiff,ωSelect = maximum(θDiv);
                    push!(divDict[i],ωSelect);
                end
            end
        end
        divSetNew,divDetNew = splitAny(divSetNew,divDetNew,Ω,divDict);
        n += 1;
    end
    return divSetNew,divDetNew;
end

function combinePart(pData,disData,Ω,divSetOld,divDetOld,H,tcoreList,xcoreList,ycoreList,wp,noTh,ubCost)
    # first combine the divSet and divDet to create new partitions
    # new partitions should be easy to solve
    # obtain the solution with the coarse partition, which will be guaranteed to be feasible
    divSet = Dict();
    divDet = Dict();
    # combine the partition
    for i in pData.II
        divSet[i] = [];
        divDet[i] = [];
        earliest0 = 0;
        latest0 = length(H) - 1;
        for par in 1:length(divSetOld[i])
            if divDetOld[i][par] == 1
                earliest0 = divSetOld[i][par].endH;
            end
            if (divDetOld[i][par] == -1)&(divSetOld[i][par].startH < latest0)
                latest0 = divSetOld[i][par].startH;
            end
        end
        # create new partition
        if earliest0 > 0
            push!(divSet[i],partType(0,earliest0));
            push!(divDet[i],1);
        end
        push!(divSet[i],partType(earliest0,latest0));
        push!(divDet[i],0);
        if latest0 < length(H) - 1
            push!(divSet[i],partType(latest0,length(H) - 1));
            push!(divDet[i],-1);
        end
    end

    ubCLTemp = [Inf];

    function partBendersC(cb)
        currentLB = MathProgBase.cbgetbestbound(cb);
        if currentLB <= minimum(ubCLTemp)
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
            # generate cuts
            θInt = Dict();
            ubTemp,θInt = ubCalP(pData,disData,Ω,xhat,that,Tmax1,1,wp);
            push!(ubCLTemp,ubTemp);

            # obtain the cores
            tcore,xcore,ycore = avgCore(pData,divSet,tcoreList,xcoreList,ycoreList);
            # here is the issue, pack it in a function prevent separating it
            dataList = subPara(pData,disData,Ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore,wp);
            cutScen = [];
            for ω in Ω
                if length(dataList[ω]) != 3
                    if (dataList[ω][4] - θhat[findfirst(Ω,ω)] > 1e-4*θhat[findfirst(Ω,ω)])
                        push!(cutScen,ω);
                    end
                end
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
        else
            return JuMP.StopTheSolver;
        end
    end

    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-8, FeasibilityTol = 1e-8, Method = 1, Threads = noTh, Cutoff = ubCost));
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

    addlazycallback(mp, partBendersC);
    mpStatus = solve(mp);
    tfeas = Dict();
    xfeas = Dict();
    for i in pData.II
        tfeas[i] = getvalue(mp[:t][i]);
        for j in pData.Ji[i]
            xfeas[i,j] = getvalue(mp[:x][i,j]);
        end
    end

    return tfeas,xfeas;
end
