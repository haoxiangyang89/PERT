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

function splitPar_CI(PartSet,PartDet,splitInfo)
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
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
                else
                    part1 = partType(PartSet[i][par].startH,splitStart);
                    part2 = partType(splitStart,PartSet[i][par].endH);
                    push!(partSetiTemp,part1);
                    push!(partSetiTemp,part2);
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
                end
            else
                if (splitEnd < PartSet[i][par].endH)&(splitEnd > PartSet[i][par].startH)
                    part1 = partType(PartSet[i][par].startH,splitEnd);
                    part2 = partType(splitEnd,PartSet[i][par].endH);
                    push!(partSetiTemp,part1);
                    push!(partSetiTemp,part2);
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
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

function iniPart(pData,disData,Ω,sN,MM)
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
    for m in 1:MM
        #randList = rand(Ω,sN);
        disData1 = Dict();
        Ω1 = 1:sN;
        randList = [(ω - 1)*MM + m for ω in Ω1];
        for i in Ω1
            disData1[i] = deepcopy(disData[randList[i]]);
            disData1[i].prDis = (1 - pData.p0)/length(Ω1);
        end
        text,xext,fext,gext,mext = extForm_cheat(pData,disData1,Ω1,1e-4,999999);
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
    end

    # return the min-max decomposition
    tHList = [];
    for i in pData.II
        tHstart = maximum([ω for ω in keys(H) if H[ω] <= minimum(tList[i])]);
        tHend = minimum([ω for ω in keys(H) if H[ω] >= maximum(tList[i])]);
        push!(tHList,[i,tHstart,tHend]);
    end
    return ubList,tHList,ubMin,tBest,xBest,θBest;
end

function splitAny(PartSet,PartDet,splitInfo)
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
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
                else
                    part1 = partType(PartSet[i][par].startH,splitStart);
                    part2 = partType(splitStart,PartSet[i][par].endH);
                    push!(partSetiTemp,part1);
                    push!(partSetiTemp,part2);
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
                end
            else
                if (splitEnd < PartSet[i][par].endH)&(splitEnd > PartSet[i][par].startH)
                    part1 = partType(PartSet[i][par].startH,splitEnd);
                    part2 = partType(splitEnd,PartSet[i][par].endH);
                    push!(partSetiTemp,part1);
                    push!(partSetiTemp,part2);
                    push!(partDetiTemp,0);
                    push!(partDetiTemp,0);
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
