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
    for (i,splitPos) in splitInfo
        partSetiTemp = [];
        partDetiTemp = [];
        for par in 1:length(PartSet[i])
            # split the current partition
            if (splitPos < PartSet[i][par].endH)&(splitPos > PartSet[i][par].startH)
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

function prunePart(pData,disData,Ω,divSet,divDet,cutSet,Tmax,currentUB)
    # try to solve the linear relaxation of the master to prune some partition
    for i in pData.II
        for par in 1:length(divSet[i])
            if divDet[i][par] == 0
                mpTemp = createMaster_DivRel(pData,disData,Ω,divSet,divDet,cutSet,Tmax,[i,par]);
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
