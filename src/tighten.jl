# This is the function to find valid inequalities given temporal relationships
# added the feature that multiple disruption could be at the same time
# added the feature that preclude some branching options

function getOmegaSeq(disData)
    sortedDis = sort(collect(disData),by = x->x[2].H);
    ωSeq = [];
    ωDict = Dict();
    currentTime = 0;
    Hω = [];
    for item in sortedDis
        if item[2].H > currentTime
            push!(ωSeq,item[1]);
            push!(Hω,disData[item[1]].H);
            ωDict[item[1]] = length(ωSeq);
        else
            if typeof(ωSeq[length(ωSeq)]) == Int64
                ωSeq[length(ωSeq)] = [ωSeq[length(ωSeq)],item[1]];
                ωDict[item[1]] = length(ωSeq);
            else
                push!(ωSeq[length(ωSeq)],item[1]);
                ωDict[item[1]] = length(ωSeq);
            end
        end
        currentTime = item[2].H;
    end
    Ωt = 1:length(Hω);
    return ωSeq,ωDict,Hω,Ωt;
end

function precludeRel(pData,disData,Ω,ub = Inf)
    # solve |II| number of LP to obtain each activity's earliest starting time
    # initialize the brInfo
    brInfo = zeros(length(pData.II),length(Ω));
    for i in pData.II
        tearly = iSolve(pData,i);
        for ω in Ω
            if tearly >= disData[ω].H
                brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] = 1;
            end
        end
    end
    if ub != Inf
        for i in pData.II
            tlate = lSolve(pData,i,ub);
            for ω in Ω
                if tlate <= disData[ω].H
                    brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] = -1;
                end
            end
        end
    end

    return brInfo;
end

function precludeRelNew(pData,H,ub = Inf)
    # solve |II| number of LP to obtain each activity's earliest starting time
    # initialize the brInfo
    brInfo = zeros(length(pData.II),length(H) - 2);
    for i in pData.II
        tearly = iSolve(pData,i);
        for hIter in keys(H)
            if (hIter != 0)&(hIter != length(H) - 1)
                if tearly >= H[hIter]
                    brInfo[findin(pData.II,i)[1],hIter] = 1;
                end
            end
        end
    end
    if ub != Inf
        for i in pData.II
            tlate = lSolve(pData,i,ub);
            for hIter in keys(H)
                if (hIter != 0)&(hIter != length(H) - 1)
                    if tlate <= H[hIter]
                        brInfo[findin(pData.II,i)[1],hIter] = -1;
                    end
                end
            end
        end
    end

    return brInfo;
end

function precludeRelStrata(pData,Hω,Ωt,ub = Inf)
    # solve |II| number of LP to obtain each activity's earliest starting time
    # initialize the brInfo
    brInfo = zeros(length(pData.II),length(Ωt));
    for i in pData.II
        tearly = iSolve(pData,i);
        for ω in Ωt
            if tearly >= Hω[ω]
                brInfo[findin(pData.II,i)[1],ω] = 1;
            end
        end
    end
    if ub != Inf
        tlate = lSolve(pData,i);
        for ω in Ωt
            if tlate <= Hω[ω]
                brInfo[findin(pData.II,i)[1],ω] = -1;
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

# exploit the logical relationships to tighten the divSet and divDet
# assume for every activity, there are only three partitions [-1,0,1]
function divExploit(pData,disData,H,PartSet,PartDet,distanceDict)
    iterBool = true;
    earlyT = Dict();
    lateT = Dict();
    for i in pData.II
        # earlyT is the index of H that is right before the earliest possible starting time of i
        # earlyT is the index of H that is right after the latest possible starting time of i
        earlyT[i] = 0;
        lateT[i] = maximum(keys(H));
        for l in 1:length(PartSet[i])
            if PartDet[i][l] == 1
                earlyT[i] = PartSet[i][l].endH;
            end
        end
        for l in length(PartSet[i]):-1:1
            if PartDet[i][l] == -1
                lateT[i] = PartSet[i][l].startH;
            end
        end
    end
    while iterBool
        # if the partition is updated then keep iteration
        iterBool = false;
        for i in pData.II
            earlyPoss = H[earlyT[i]];
            latePoss = H[lateT[i]];
            for j in pData.II
                if (i,j) in keys(distanceDict)
                    # if j is after i
                    if latePoss > H[maximum([PartSet[j][l].endH for l in 1:length(PartSet[j]) if PartDet[j][l] == 0])] - distanceDict[i,j]
                        latePoss = H[maximum([PartSet[j][l].endH for l in 1:length(PartSet[j]) if PartDet[j][l] == 0])] - distanceDict[i,j];
                    end
                elseif (j,i) in keys(distanceDict)
                    # if j is before i
                    if earlyPoss < H[minimum([PartSet[j][l].startH for l in 1:length(PartSet[j]) if PartDet[j][l] == 0])] + distanceDict[j,i]
                        earlyPoss = H[minimum([PartSet[j][l].startH for l in 1:length(PartSet[j]) if PartDet[j][l] == 0])] + distanceDict[j,i];
                    end
                end
            end
            lateSet = [hInd for hInd in keys(H) if H[hInd] > latePoss];
            if lateSet != []
                lateTNew = minimum(lateSet);
            else
                lateTNew = maximum(keys(H));
            end
            earlySet = [hInd for hInd in keys(H) if H[hInd] <= earlyPoss];
            earlyTNew = maximum(earlySet);
            if lateTNew < lateT[i]
                lateT[i] = lateTNew;
                iterBool = true;
            end
            if earlyTNew > earlyT[i]
                earlyT[i] = earlyTNew;
                iterBool = true;
            end
        end
    end

    divSetNew = Dict();
    divDetNew = Dict();
    for i in pData.II
        divSetNew[i] = [];
        divDetNew[i] = [];
        for l in 1:length(PartSet[i])
            # 1 set
            if earlyT[i] >= PartSet[i][l].endH
                push!(divSetNew[i],deepcopy(PartSet[i][l]));
                push!(divDetNew[i],1);
            elseif (earlyT[i] > PartSet[i][l].startH)&(earlyT[i] < PartSet[i][l].endH)&(PartDet[i][l] == 0)
                divSetNew[i][length(divSetNew[i])].endH = earlyT[i];
                push!(divSetNew[i],deepcopy(PartSet[i][l]));
                divSetNew[i][length(divSetNew[i])].startH = earlyT[i];
                push!(divDetNew[i],0);
            else
                push!(divSetNew[i],deepcopy(PartSet[i][l]));
                push!(divDetNew[i],0);
            end
        end
        for l in length(divSetNew[i]):-1:1
            # -1 set
            if lateT[i] <= divSetNew[i][l].startH
                divDetNew[i][l] = -1;
            elseif (lateT[i] >= divSetNew[i][l].startH)&(lateT[i] < divSetNew[i][l].endH)&(divDetNew[i][l] == 0)
                if l + 1 <= length(divSetNew[i])
                    # if there is another partition after this
                    divSetNew[i][l+1].startH = lateT[i];
                    divSetNew[i][l].endH = lateT[i];
                else
                    # if this is the last partition
                    push!(divSetNew[i],deepcopy(PartSet[i][l]));
                    divSetNew[i][l].endH = lateT[i];
                    divSetNew[i][l+1].startH = lateT[i];
                    push!(divDetNew[i],-1);
                end
            end
        end
    end
    return divSetNew,divDetNew;
end

# break the divSet and get 2 updated divSets and divDets
function breakDiv(pData,disData,H,divSet,divDet,iBreak,locBreak,distanceDict)
    # iBreak is the activity to branch on
    # locBreak is the point of H to break on

    # obtain the branching pieces
    stopBool = false;
    currentDiv = 1;
    while !(stopBool)
        if (locBreak < divSet[iBreak][currentDiv].endH)&(locBreak >= divSet[iBreak][currentDiv].startH)
            stopBool = true;
        else
            currentDiv += 1;
        end
    end

    newPartSet1 = deepcopy(divSet);
    newPartDet1 = deepcopy(divDet);
    newPartSet1[iBreak] = Array{partType,1}();
    newPartDet1[iBreak] = Array{Int64,1}();
    for l in 1:length(divSet[iBreak])
        if l == currentDiv
            push!(newPartSet1[iBreak],partType(divSet[iBreak][l].startH,locBreak));
            push!(newPartDet1[iBreak],1);
            push!(newPartSet1[iBreak],partType(locBreak,divSet[iBreak][l].endH));
            push!(newPartDet1[iBreak],0);
        else
            push!(newPartSet1[iBreak],divSet[iBreak][l]);
            push!(newPartDet1[iBreak],divDet[iBreak][l]);
        end
    end
    divSet1,divDet1 = divExploit(pData,disData,H,newPartSet1,newPartDet1,distanceDict);

    newPartSet2 = deepcopy(divSet);
    newPartDet2 = deepcopy(divDet);
    newPartSet2[iBreak] = Array{partType,1}();
    newPartDet2[iBreak] = Array{Int64,1}();
    for l in 1:length(divSet[iBreak])
        if l == currentDiv
            push!(newPartSet2[iBreak],partType(divSet[iBreak][l].startH,locBreak));
            push!(newPartDet2[iBreak],0);
            push!(newPartSet2[iBreak],partType(locBreak,divSet[iBreak][l].endH));
            push!(newPartDet2[iBreak],-1);
        else
            push!(newPartSet2[iBreak],divSet[iBreak][l]);
            push!(newPartDet2[iBreak],divDet[iBreak][l]);
        end
    end
    divSet2,divDet2 = divExploit(pData,disData,H,newPartSet2,newPartDet2,distanceDict);

    return divSet1,divDet1,divSet2,divDet2;
end

function obtainBds(pData,disData,Ω,mpTemp,ub)
    H = Dict();
    H[0] = 0;
    for ω in Ω
        H[ω] = disData[ω].H;
    end
    # obtain the boundInfo
    ubInfo = Dict();
    lbInfo = Dict();
    for i in pData.II
        ubInfo[i] = disData[length(Ω)].H + sum(pData.D[i] for i in pData.II);
        lbInfo[i] = 0;
    end

    keepIter = true;
    @constraint(mpTemp,pData.p0*mpTemp[:t][0] + sum(disData[ω].prDis*mpTemp[:θ][ω] for ω in Ω) <= ub);
    ubTempInfo = Dict();
    lbTempInfo = Dict();
    while keepIter
        keepIter = false;
        for i in pData.II
            @constraint(mpTemp,mpTemp[:t][i] <= ubInfo[i]);
            @constraint(mpTemp,mpTemp[:t][i] >= lbInfo[i]);
            @objective(mpTemp,Max,mpTemp[:t][i]);
            solve(mpTemp);
            ubTempInfo[i] = getobjectivevalue(mpTemp);
            @objective(mpTemp,Min,mpTemp[:t][i]);
            solve(mpTemp);
            lbTempInfo[i] = getobjectivevalue(mpTemp);
            if (ubTempInfo[i] < ubInfo[i] - 1e-5)
                ubInfo[i] = ubTempInfo[i];
                keepIter = true;
            end
            if (lbTempInfo[i] > lbInfo[i] + 1e-5)
                lbInfo[i] = lbTempInfo[i];
                keepIter = true;
            end
        end
    end
    return ubInfo,lbInfo;
end

function longestPath(pData,pInfo = Dict(),outLink = false)
    lDict = Dict();
    finishedList = [];
    activeList = [];
    linkDict = Dict();
    for i in pData.II
        if pData.Pre[i] == []
            lDict[i] = 0;
        else
            lDict[i] = -Inf;
        end
        linkDict[i] = -1;
    end
    if pInfo == Dict()
        pd = pData.D;
    else
        pd = pInfo;
    end
    # Dijkstra with negative edge weight
    while length(finishedList) < length(pData.II)
        for i in pData.II
            if !((i in activeList)|(i in finishedList))
                enterBool = true;
                for j in pData.Pre[i]
                    if !(j in finishedList)
                        enterBool = false;
                    end
                end
                if enterBool
                    push!(activeList,i);
                end
            end
        end
        # find the largest activity in activeList
        maxVal = -Inf;
        maxInd = -1;
        for iInd in activeList
            if maxVal < lDict[iInd]
                maxVal = lDict[iInd];
                maxInd = iInd;
            end
        end
        # update the activities connected to maxInd
        for j in pData.Succ[maxInd]
            if lDict[j] < lDict[maxInd] + pd[maxInd]
                lDict[j] = lDict[maxInd] + pd[maxInd];
                linkDict[j] = maxInd;
            end
        end
        push!(finishedList,maxInd);
        deleteat!(activeList,findfirst(activeList,maxInd));
    end
    if outLink
        return lDict,linkDict;
    else
        return lDict;
    end
end

function obtainDet(pData,disData,Ω,mpTemp,ub,divSet,divDet)
    # predetermine some y's data
    @constraint(mpTemp,pData.p0*mpTemp[:t][0] + sum(disData[ω].prDis*mpTemp[:θ][ω] for ω in Ω) <= ub);
    for i in pData.II
        # fix the y[i,par] to be 1 and solve the problem, if infeasible, then y[i,par] has to be 0
        for par in 1:length(divSet[i])
            if divDet[i][par] == 0
                mpTTemp = copy(mpTemp);
                @constraint(mpTTemp,mpTTemp[:y][i,par] == 1);
                mpStatus = solve(mpTTemp);
                if mpStatus == :Infeasible
                    divDet[i][par] = -1;
                end
            end
        end
    end
    return divDet;
end

function detCal(pData,ii,jj)
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variables(mp, begin
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
    end);
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @objective(mp, Min, t[jj] - t[ii]);
    solve(mp);
    disIJ = getobjectivevalue(mp);

    return disIJ;
end

function closestCore(pData,divSet,ycoreTemp,tcoreList,xcoreList,ycoreList)
    smallestDist = Inf;
    minInd = -1;
    for ll in 1:length(ycoreList)
        distanceCore = sum(sum(abs(ycoreList[ll][i,par] - ycoreTemp[i,par]) for par in 1:length(divSet[i])) for i in pData.II);
        if distanceCore < smallestDist
            minInd = ll;
            smallestDist = distanceCore;
        end
    end
    return tcoreList[minInd],xcoreList[minInd],ycoreList[minInd];
end

function avgCore(pData,divSet,tcoreList,xcoreList,ycoreList)
    tcore = Dict();
    xcore = Dict();
    ycore = Dict();
    for i in pData.II
        tcore[i] = mean([tcoreList[ll][i] for ll in 1:length(tcoreList)]);
        for j in pData.Ji[i]
            xcore[i,j] = mean([xcoreList[ll][i,j] for ll in 1:length(xcoreList)]);
        end
        for par in 1:length(divSet[i])
            ycore[i,par] = mean([ycoreList[ll][i,par] for ll in 1:length(ycoreList)]);
        end
    end
    return tcore,xcore,ycore;
end

function avgCoreShare(pData,H,divSet,tcoreP,xcoreP,tcoreN,xcoreN,weightP)
    tcoreO = Dict();
    xcoreO = Dict();
    ycoreO = Dict();
    for i in pData.II
        tcoreO[i] = (tcoreP[i]*weightP+tcoreN[i])/(weightP+1);
        for j in pData.Ji[i]
            xcoreO[i,j] = (xcoreP[i,j]*weightP+xcoreN[i,j])/(weightP+1);
        end
        for par in 1:length(divSet[i])
            if (tcoreO[i] >= H[divSet[i][par].startH])&(tcoreO[i] < H[divSet[i][par].endH])
                ycoreO[i,par] = 1;
            elseif (abs(tcoreO[i] - H[length(H) - 1]) < 1e-4)&(divSet[i][par].endH == length(H) - 1)
                ycoreO[i,par] = 1;
            else
                ycoreO[i,par] = 0;
            end
        end
    end
    return tcoreO,xcoreO,ycoreO;
end

function genCritical(pData,lInfo)
    lGGDict,linkGDict = longestPath(pData,lInfo,true);
    currentA = 0;
    criticalList = [];
    while currentA != -1
        push!(criticalList,currentA);
        currentA = linkGDict[currentA];
    end
    return criticalList;
end
