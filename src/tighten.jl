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

function longestPath(pData)
    lDict = Dict();
    finishedList = [];
    activeList = [];
    for i in pData.II
        if pData.Pre[i] == []
            lDict[i] = 0;
        else
            lDict[i] = -Inf;
        end
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
            if lDict[j] < lDict[maxInd] + pData.D[maxInd]
                lDict[j] = lDict[maxInd] + pData.D[maxInd];
            end
        end
        push!(finishedList,maxInd);
        deleteat!(activeList,findfirst(activeList,maxInd));
    end
    return lDict;
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

function detCal(pData,i,j)
    mp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variables(mp, begin
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
    end);
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    @objective(mp, Min, t[j] - t[i]);
    solve(mp);
    disIJ = getobjectivevalue(mp);

    return disIJ;
end
