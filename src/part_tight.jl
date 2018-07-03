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
    newPartSet = copy(PartSet);
    newPartDet = copy(PartDet);
    # for each activity
    for i in pData.II
        partSetiTemp = [];
        partDetiTemp = [];
        for par in 1:length(PartSet[i])
            if PartDet[i][par] == 0
                if (PartSet[i][par].startH < lbInfo[i])
                    # split the current set into parts with the first det as 1
                    set1start = PartSet[i][par].startH;
                    set1end = maximum([ω for ω in PartSet[i][par].startH:PartSet[i][par].endH if disData[ω].H < lbInfo[i]]);
                    if set1end == PartSet[i][par].endH
                        push!(partSetiTemp,PartSet[i][par]);
                        push!(partDetiTemp,1);
                    else
                        push!(partSetiTemp,partType(set1start,set1end));
                        push!(partDetiTemp,1);
                        set2start = set1end;
                        if (PartSet[i][par].endH < ubInfo[i])
                            set2end = PartSet[i][par].endH;
                            push!(partSetiTemp,partType(set2start,set2end));
                            push!(partDetiTemp,0);
                        else
                            set2end = minimum([ω for ω in PartSet[i][par].startH:PartSet[i][par].endH if disData[ω].H >= ubInfo[i]]);
                            push!(partSetiTemp,partType(set2start,set2end));
                            push!(partDetiTemp,0);
                            set3start = set2end;
                            set3end = PartSet[i][par].endH;
                            push!(partSetiTemp,partType(set3start,set3end));
                            push!(partDetiTemp,-1);
                        end
                    end
                    # split the current set into parts with the first det as 0
                elseif (PartSet[i][par].startH < ubInfo[i])
                    set1start = PartSet[i][par].startH;
                    set1end = minimum([ω for ω in PartSet[i][par].startH:PartSet[i][par].endH if disData[ω].H >= ubInfo[i]]);
                    if set1end == PartSet[i][par].endH
                        push!(partSetiTemp,PartSet[i][par]);
                        push!(partDetiTemp,0);
                    else
                        set2start = set1end;
                        set2end = PartSet[i][par].endH;
                        push!(partSetiTemp,partType(set2start,set2end));
                        push!(partDetiTemp,-1);
                    end
                else
                    push!(partSetiTemp,PartSet[i][par]);
                    push!(partDetiTemp,-1);
                end
            else
                push!(partSetiTemp,PartSet[i][par]);
                push!(partDetiTemp,PartDet[i][par]);
            end
        end
        newPartSet[i] = partSetiTemp;
        newPartDet[i] = partDetiTemp;
    end
end
