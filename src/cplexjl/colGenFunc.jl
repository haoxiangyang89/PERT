# create a master problem with a selection of binary variables
function masterMax(pData,disData,Ω,ωInfoOld,Ghat,ωInfoNew,cutSet,partCurrentTemp,partDetTemp,M = 9999999, solveOpt = 0)
    partRev = Dict();
    partNo = Dict();
    for i in pData.II
        partRev[i] = Dict();
        partNo[i] = length(partCurrent[i]);
        for partIter in 1:partNo[i]
            for item in partCurrent[i][partIter]
                partRev[i][item] = partIter;
            end
        end
    end
    # cutSet includes the cuts generated
    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-9,FeasibilityTol = 1e-9,OutputFlag = 0));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
      0 <= G[i in pData.II,ω in Ω] <= 1
    end);
    for (i,ω) in ωInfoNew
        setcategory(G[i,ω], :Bin);
    end
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    # logic constraints between G and t
    @constraint(mp, tGcons1I[i in pData.II, ω in Ω], t[i] <= disData[ω].H + G[i,ω]*M);
    @constraint(mp, tGcons2I[i in pData.II, ω in Ω], t[i] >= disData[ω].H - (1 - G[i,ω])*M);
    # logic constraints of G
    @constraint(mp, GijCon[i in pData.II, j in pData.Succ[i],ω in Ω], G[i,ω] <= G[j,ω]);
    @constraint(mp, Gω12Con[i in pData.II, ω in 1:(length(Ω)-1)], G[i,ω] >= G[i,ω + 1]);
    @constraint(mp, Gdet[(i,ω) in ωInfoOld], G[i,ω] == Ghat[i,ω]);
    @constraint(mp, Gdet1[i in pData.II, ω in Ω; partDet[i][partRev[i][ω]] == 1], G[i,ω] == 1);
    @constraint(mp, Gdet0[i in pData.II, ω in Ω; partDet[i][partRev[i][ω]] == -1], G[i,ω] == 0);

    for ω in Ω
        if cutSet[ω] != []
            for nc in 1:length(cutSet[ω])
                @constraint(mp, θ[ω] >= cutSet[ω][nc][4] + sum(cutSet[ω][nc][1][i]*(mp[:t][i] - cutSet[ω][nc][5][i]) for i in pData.II) +
                    sum(sum(cutSet[ω][nc][2][i,j]*(mp[:x][i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II) +
                    sum(cutSet[ω][nc][3][i]*(mp[:G][i,ω] - cutSet[ω][nc][7][i]) for i in pData.II));
            end
        end
    end

    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));
    if solveOpt == 1
        solve(mp);
        return getobjectivevalue(mp);
    else
        return mp;
    end
end

function masterMaxLP(pData,disData,Ω,ωInfoOld,Ghat,ωInfoNew,cutSet,partCurrentTemp,partDetTemp,M = 9999999)
    partRev = Dict();
    partNo = Dict();
    for i in pData.II
        partRev[i] = Dict();
        partNo[i] = length(partCurrent[i]);
        for partIter in 1:partNo[i]
            for item in partCurrent[i][partIter]
                partRev[i][item] = partIter;
            end
        end
    end
    # cutSet includes the cuts generated
    mp = Model(solver = GurobiSolver(IntFeasTol = 1e-9,FeasibilityTol = 1e-9,OutputFlag = 0));
    @variables(mp, begin
      θ[Ω] >= 0
      0 <= x[i in pData.II,j in pData.Ji[i]] <= 1
      t[i in pData.II] >= 0
      0 <= G[i in pData.II,ω in Ω] <= 1
    end);
    @constraint(mp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*x[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);

    # logic constraints between G and t
    @constraint(mp, tGcons1I[i in pData.II, ω in Ω], t[i] <= disData[ω].H + G[i,ω]*M);
    @constraint(mp, tGcons2I[i in pData.II, ω in Ω], t[i] >= disData[ω].H - (1 - G[i,ω])*M);
    # logic constraints of G
    @constraint(mp, GijCon[i in pData.II, j in pData.Succ[i],ω in Ω], G[i,ω] <= G[j,ω]);
    @constraint(mp, Gω12Con[i in pData.II, ω in 1:(length(Ω)-1)], G[i,ω] >= G[i,ω + 1]);
    @constraint(mp, Gdet[(i,ω) in ωInfoOld], G[i,ω] == Ghat[i,ω]);
    @constraint(mp, Gdet1[i in pData.II, ω in Ω; partDet[i][partRev[i][ω]] == 1], G[i,ω] == 1);
    @constraint(mp, Gdet0[i in pData.II, ω in Ω; partDet[i][partRev[i][ω]] == -1], G[i,ω] == 0);

    for ω in Ω
        if cutSet[ω] != []
            for nc in 1:length(cutSet[ω])
                @constraint(mp, θ[ω] >= cutSet[ω][nc][4] + sum(cutSet[ω][nc][1][i]*(mp[:t][i] - cutSet[ω][nc][5][i]) for i in pData.II) +
                    sum(sum(cutSet[ω][nc][2][i,j]*(mp[:x][i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II) +
                    sum(cutSet[ω][nc][3][i]*(mp[:G][i,ω] - cutSet[ω][nc][7][i]) for i in pData.II));
            end
        end
    end

    @objective(mp, Min, pData.p0*t[0] + sum(disData[ω].prDis*θ[ω] for ω in Ω));
    for (i,ω) in ωInfoNew
        JuMP.fix(G[i,ω],0);
    end
    solve(mp);
    obj0 = getobjectivevalue(mp);
    for (i,ω) in ωInfoNew
        JuMP.fix(G[i,ω],1);
    end
    solve(mp);
    obj1 = getobjectivevalue(mp);
    if obj0 < obj1
        finalobj = obj0;
    else
        finalobj = obj1;
    end
    return finalobj;
end

function searchBest(pData,disData,Ω,ωInfo,Ghat,cutSet,partCurrentTemp,partDetTemp,M = 999999)
    # start with activities with no predecessors
    activityList = [];
    finishedList = [];
    while length(finishedList) != length(pData.II)
        # obtain the activity to search for
        for i in pData.II
            if (!(i in finishedList))&(!(i in activityList))
                addBool = true;
                for j in pData.Pre[i]
                    if !(j in finishedList)
                        addBool = false;
                    end
                end
                if addBool
                    push!(activityList,i);
                end
            end
        end
        currenti = activityList[1];
        shift!(activityList);

        # search process
        fracList = [ω for ω in Ω if (Ghat[currenti,ω] > 1e-7)&(Ghat[currenti,ω] < 1 - 1e-7)];
        if fracList != []
            startInd = minimum(fracList);
            endInd = maximum(fracList);
            halfInd = Int64(ceil((startInd + endInd)/2));
            startObj = masterMaxLP(pData,disData,Ω,ωInfo,Ghat,[(currenti,startInd)],cutSet,partCurrentTemp,partDetTemp,M);
            endObj = masterMaxLP(pData,disData,Ω,ωInfo,Ghat,[(currenti,endInd)],cutSet,partCurrentTemp,partDetTemp,M);
            halfObj = masterMaxLP(pData,disData,Ω,ωInfo,Ghat,[(currenti,halfInd)],cutSet,partCurrentTemp,partDetTemp,M);
            while endInd - startInd > 2
                q1Ind = Int64(ceil((startInd + halfInd)/2));
                q3Ind = Int64(ceil((halfInd + endInd)/2));
                q1Obj = masterMaxLP(pData,disData,Ω,ωInfo,Ghat,[(currenti,q1Ind)],cutSet,partCurrentTemp,partDetTemp,M);
                q3Obj = masterMaxLP(pData,disData,Ω,ωInfo,Ghat,[(currenti,q3Ind)],cutSet,partCurrentTemp,partDetTemp,M);
                if (round(startObj,7) > round(q1Obj,7))|(round(q1Obj,7) > round(halfObj))
                    halfIndTemp = halfInd;
                    halfObjTemp = halfObj;
                    halfInd = q1Ind;
                    halfObj = q1Obj;
                    endInd = halfIndTemp;
                    endObj = halfObjTemp;
                elseif (round(q1Obj,7) > round(halfObj,7))|(round(halfObj,7) > round(q3Obj,7))
                    startInd = q1Ind;
                    startObj = q1Obj;
                    endInd = q3Ind;
                    endObj = q3Obj;
                else
                    halfIndTemp = halfInd;
                    halfObjTemp = halfObj;
                    halfInd = q3Ind;
                    halfObj = q3Obj;
                    startInd = halfIndTemp;
                    startObj = halfObjTemp;
                end
            end
            if (startObj >= halfObj)&(halfObj >= endObj)
                # startInd is the best
                mpstart = masterMax(pData,disData,Ω,ωInfo,Ghat,[(currenti,startInd)],cutSet,partCurrentTemp,partDetTemp,M);
                solve(mpstart);
                GhatTemp = getvalue(mpstart[:G]);
                for i in pData.II
                    for ω in Ω
                        if GhatTemp[i,ω] < 1e-7
                            Ghat[i,ω] = 0;
                        elseif GhatTemp[i,ω] > 1 - 1e-7
                            Ghat[i,ω] = 1;
                        else
                            Ghat[i,ω] = GhatTemp[i,ω];
                        end
                    end
                end
                push!(ωInfo,(currenti,startInd));
            elseif (startObj <= halfObj)&(halfObj >= endObj)
                # halfInd is the best
                mphalf = masterMax(pData,disData,Ω,ωInfo,Ghat,[(currenti,halfInd)],cutSet,partCurrentTemp,partDetTemp,M);
                solve(mphalf);
                GhatTemp = getvalue(mphalf[:G]);
                for i in pData.II
                    for ω in Ω
                        if GhatTemp[i,ω] < 1e-7
                            Ghat[i,ω] = 0;
                        elseif GhatTemp[i,ω] > 1 - 1e-7
                            Ghat[i,ω] = 1;
                        else
                            Ghat[i,ω] = GhatTemp[i,ω];
                        end
                    end
                end
                push!(ωInfo,(currenti,halfInd));
            else
                # halfInd is the best
                mpend = masterMax(pData,disData,Ω,ωInfo,Ghat,[(currenti,endInd)],cutSet,partCurrentTemp,partDetTemp,M);
                solve(mpend);
                GhatTemp = getvalue(mpend[:G]);
                for i in pData.II
                    for ω in Ω
                        if GhatTemp[i,ω] < 1e-7
                            Ghat[i,ω] = 0;
                        elseif GhatTemp[i,ω] > 1 - 1e-7
                            Ghat[i,ω] = 1;
                        else
                            Ghat[i,ω] = GhatTemp[i,ω];
                        end
                    end
                end
                push!(ωInfo,(currenti,endInd));
            end
        end
        push!(finishedList,currenti);
    end
    return ωInfo;
end

function searchBestnou(pData,disData,Ω,ωInfo,Ghat,cutSet,partCurrentTemp,partDetTemp,M = 999999)
    # start with activities with no predecessors
    activityList = [];
    finishedList = [];
    ωInfoNew = [];
    while length(finishedList) != length(pData.II)
        # obtain the activity to search for
        for i in pData.II
            if (!(i in finishedList))&(!(i in activityList))
                addBool = true;
                for j in pData.Pre[i]
                    if !(j in finishedList)
                        addBool = false;
                    end
                end
                if addBool
                    push!(activityList,i);
                end
            end
        end
        currenti = activityList[1];
        shift!(activityList);

        # search process
        fracList = [ω for ω in Ω if (Ghat[currenti,ω] > 1e-7)&(Ghat[currenti,ω] < 1 - 1e-7)];
        if fracList != []
            startInd = minimum(fracList);
            endInd = maximum(fracList);
            halfInd = Int64(ceil((startInd + endInd)/2));
            startObj = masterMaxLP(pData,disData,Ω,ωInfo,Ghat,[(currenti,startInd)],cutSet,partCurrentTemp,partDetTemp,M);
            endObj = masterMaxLP(pData,disData,Ω,ωInfo,Ghat,[(currenti,endInd)],cutSet,partCurrentTemp,partDetTemp,M);
            halfObj = masterMaxLP(pData,disData,Ω,ωInfo,Ghat,[(currenti,halfInd)],cutSet,partCurrentTemp,partDetTemp,M);
            while endInd - startInd > 2
                q1Ind = Int64(ceil((startInd + halfInd)/2));
                q3Ind = Int64(ceil((halfInd + endInd)/2));
                q1Obj = masterMaxLP(pData,disData,Ω,ωInfo,Ghat,[(currenti,q1Ind)],cutSet,partCurrentTemp,partDetTemp,M);
                q3Obj = masterMaxLP(pData,disData,Ω,ωInfo,Ghat,[(currenti,q3Ind)],cutSet,partCurrentTemp,partDetTemp,M);
                objList = [startObj,q1Obj,halfObj,q3Obj,endObj];
                # one of 15 possibilities
                if round(startObj,7) == round(maximum(objList),7)
                    if round(q3Obj,7) < round(startObj,7)
                        scenNo = 1;
                    else
                        scenNo = 2;
                    end
                elseif round(q1Obj,7) == round(maximum(objList),7)
                    if round(halfObj,7) < round(q1Obj,7)
                        scenNo = 1;
                    else
                        scenNo = 2;
                    end
                elseif round(halfObj,7) == round(maximum(objList),7)
                    if round(endObj,7) < round(halfObj,7)
                        scenNo = 2;
                    else
                        scenNo = 3;
                    end
                else
                    scenNo = 3;
                end
                # ((round(startObj,7) == round(q1Obj,7))&(round(q1Obj,7) == round(halfObj,7)))
                # ((round(q1Obj,7) == round(halfObj,7))&(round(halfObj,7) == round(q3Obj,7)))
                if scenNo == 1
                    halfIndTemp = halfInd;
                    halfObjTemp = halfObj;
                    halfInd = q1Ind;
                    halfObj = q1Obj;
                    endInd = halfIndTemp;
                    endObj = halfObjTemp;
                elseif scenNo == 2
                    startInd = q1Ind;
                    startObj = q1Obj;
                    endInd = q3Ind;
                    endObj = q3Obj;
                else
                    halfIndTemp = halfInd;
                    halfObjTemp = halfObj;
                    halfInd = q3Ind;
                    halfObj = q3Obj;
                    startInd = halfIndTemp;
                    startObj = halfObjTemp;
                end
            end
            if (round(startObj,7) >= round(halfObj,7))&(round(halfObj,7) >= round(endObj,7))
                # startInd is the best
                push!(ωInfoNew,(currenti,startInd));
            elseif (round(startObj,7) <= round(halfObj,7))&(round(halfObj,7) >= round(endObj,7))
                # halfInd is the best
                push!(ωInfoNew,(currenti,halfInd));
            else
                # endInd is the best
                push!(ωInfoNew,(currenti,endInd));
            end
        end
        push!(finishedList,currenti);
    end
    return ωInfoNew;
end
