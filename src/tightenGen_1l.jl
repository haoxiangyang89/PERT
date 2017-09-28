# This is the function to find valid inequalities given temporal relationships
# added the feature that multiple disruption could be at the same time
# added the feature that preclude some branching options

function getOmegaSeq(disData)
    sortedDis = sort(collect(disData),by = x->x[2].H);
    ωSeq = [];
    currentTime = 0;
    for item in sortedDis
        if item[2].H > currentTime
            push!(ωSeq,item[1]);
        else
            if typeof(ωSeq[length(ωSeq)]) == Int64
                ωSeq[length(ωSeq)] = [ωSeq[length(ωSeq)],item[1]];
            else
                push!(ωSeq[length(ωSeq)],item[1]);
            end
        end
        currentTime = item[2].H;
    end
    return ωSeq;
end

function precludeRel(pData,disData,Ω)
    # solve |II| number of LP to obtain each activity's earliest starting time
    # initialize the brInfo
    brInfo = zeros(length(pData.II),length(Ω));
    for i in pData.II
        tearly = iSolve(pData,disData,i);
        for ω in Ω
            if tearly >= disData[ω].H
                brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] = 1;
            end
        end
    end

    return brInfo;
end

# detect the inconsistency within the same activity
function findWrongI(i,brInfo,ωSeq)
    # no -1 should be on the left hand side of 1
    correctBool = true;
    oneRight = -1;
    negoneLeft = length(ωSeq) + 1;
    for j in 1:length(ωSeq)
        if typeof(ωSeq[j]) == Int64
            if brInfo[findin(pData.II,i)[1],findin(Ω,ωSeq[j])[1]] == -1
                negoneLeft = min(j,negoneLeft);
            end
            if brInfo[findin(pData.II,i)[1],findin(Ω,ωSeq[j])[1]] == 1
                oneRight = max(j,oneRight);
            end
        else
            for jitem in ωSeq[j]
                if brInfo[findin(pData.II,i)[1],findin(Ω,jitem)[1]] == -1
                    negoneLeft = min(j,negoneLeft);
                end
                if brInfo[findin(pData.II,i)[1],findin(Ω,jitem)[1]] == 1
                    oneRight = max(j,oneRight);
                end
            end
        end
    end
    if oneRight > negoneLeft
        correctBool = false;
    end
    return correctBool,oneRight,negoneLeft;
end

# exploint the brInfo
function brInfoExt(pData,disData,Ω,brInfo,ωSeq)
    # currently just guarantee the consistency in:
    #   the same activity: if activity i starts after ω1, then it should also start after any disruption time before disData[ω1].H
    #   the adjacent activity: if activity i starts before ω1, then all its predecessors should start before ω1
    oneRight = Dict();
    negoneLeft = Dict();
    inconBool = true;
    for i in pData.II
        correctBool,oneRight[i],negoneLeft[i] = findWrongI(i,brInfo,ωSeq);
        if correctBool
            # if there is no inconsistency, then change the brInfo
            for j in 1:oneRight[i]
                if typeof(ωSeq[j]) == Int64
                    brInfo[findin(pData.II,i)[1],findin(Ω,ωSeq[j])[1]] = 1;
                else
                    for itemj in ωSeq[j]
                        brInfo[findin(pData.II,i)[1],findin(Ω,itemj)[1]] = 1;
                    end
                end
            end
            for j in negoneLeft[i]:length(ωSeq)
                if typeof(ωSeq[j]) == Int64
                    brInfo[findin(pData.II,i)[1],findin(Ω,ωSeq[j])[1]] = -1;
                else
                    for itemj in ωSeq[j]
                        brInfo[findin(pData.II,i)[1],findin(Ω,itemj)[1]] = -1;
                    end
                end
            end
        else
            # if there is an inconsistency, then return [] for brInfo
            brInfo = [];
            inconBool = false;
            break;
        end
    end
    # test the inconsistency among precedence relationships
    if inconBool
        brInfoPrev = [];
        while brInfoPrev != brInfo
            brInfoPrev = deepcopy(brInfo);
            for i in pData.II
                # its predecessors oneRight have to be smaller than its negoneLeft
                for k in pData.Pre[i]
                    if oneRight[k] >= negoneLeft[i]
                        brInfo = [];
                        inconBool = false;
                        break;
                    end
                    # for the same scenario, the predecessors have the same branching option if the current option is -1
                    for j in Ω
                        if brInfo[findin(pData.II,i)[1],findin(Ω,j)[1]] == -1
                            if brInfo[findin(pData.II,k)[1],findin(Ω,j)[1]] <= 0
                                brInfo[findin(pData.II,k)[1],findin(Ω,j)[1]] = -1;
                            else
                                brInfo = [];
                                inconBool = false;
                                break;
                            end
                        end
                    end
                end
                # its successors negoneLeft have to be larger than its oneRight
                for k in pData.Succ[i]
                    if negoneLeft[k] <= oneRight[i]
                        brInfo = [];
                        inconBool = false;
                        break;
                    end
                    # for the same scenario, the successors have the same branching option if the current option is 1
                    for j in Ω
                        if brInfo[findin(pData.II,i)[1],findin(Ω,j)[1]] == 1
                            if brInfo[findin(pData.II,k)[1],findin(Ω,j)[1]] >= 0
                                brInfo[findin(pData.II,k)[1],findin(Ω,j)[1]] = 1;
                            else
                                brInfo = [];
                                inconBool = false;
                                break;
                            end
                        end
                    end
                end
            end
            # for the scenarios bundled together make them consistent
            for i in pData.II
                for j in 1:length(ωSeq)
                    # if it is a bundle
                    if typeof(ωSeq[j]) != Int64
                        posBool = false;
                        negBool = false;
                        for jitem in j
                            if brInfo[findin(pData.II,i)[1],findin(Ω,jitem)[1]] == 1
                                posBool = true;
                            end
                            if brInfo[findin(pData.II,i)[1],findin(Ω,jitem)[1]] == -1
                                negBool = true;
                            end
                        end
                        if (posBool)&(negBool)
                            brInfo = [];
                            inconBool = false;
                            break;
                        else
                            if posBool
                                for jitem in j
                                    brInfo[findin(pData.II,i)[1],findin(Ω,jitem)[1]] = 1;
                                end
                            end
                            if negBool
                                for jitem in j
                                    brInfo[findin(pData.II,i)[1],findin(Ω,jitem)[1]] = -1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return brInfo;
end
