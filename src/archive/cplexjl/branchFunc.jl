# this is the branching execution of the procedure
function branchSimple(pData,disData,Ω,node,that,ωSeq)
    # simply branch the smallest activity-disruption time pair
    minDist = Inf;
    minI = -1;
    minω = -1;
    for i in pData.II
        for ω in Ω
            # only branch on the non determined activity-disruption pair
            if node.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                distHt = abs(disData[ω].H - that[i]);
                if distHt < minDist
                    minDist = distHt;
                    minI = i;
                    minω = ω;
                end
            end
        end
    end
    # if we find an activity-disruption pair
    if (minI != -1)&(minω != -1)
        brInfo1 = copy(node.brInfo);
        brInfo1[findin(pData.II,minI)[1],findin(Ω,minω)[1]] = -1;
        brInfo1 = brInfoExt(pData,disData,Ω,brInfo1,ωSeq);
        brInfo2 = copy(node.brInfo);
        brInfo2[findin(pData.II,minI)[1],findin(Ω,minω)[1]] = 1;
        brInfo2 = brInfoExt(pData,disData,Ω,brInfo2,ωSeq);
    end
    # add the constraints of node1 and node2
    mp1 = copy(node.mp);
    if brInfo1 != []
        # inherit the lbCost of the branched node
        node1 = nodeType(node.lbCost,mp1,brInfo1);
    end
    mp2 = copy(node.mp);
    if brInfo2 != []
        # inherit the lbCost of the branched node
        node2 = nodeType(node.lbCost,mp2,brInfo2);
    end
    return node1,node2;
end

function branchCount(pData,disData,Ω,node,that,ωSeq)
    # branch the activity which has the most 0 in its row of brInfo
    noZeros = Dict();
    iMax = -1;
    maxZero = 0;
    for i in pData.II
        noZeros[i] = sum([node.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0 for ω in Ω]);
        if noZeros[i] > maxZero
            maxZero = noZeros[i];
            iMax = i;
        end
    end

    # find the mid point of scenarios we branch on
    first0 = 0;
    last0 = 0;
    for ω in ωSeq
       if nCurrent.brInfo[findin(pData.II,iMax)[1],findin(Ω,ω)[1]] > 0
           first0 += 1;
       end
       if nCurrent.brInfo[findin(pData.II,iMax)[1],findin(Ω,ω)[1]] > -1
           last0 += 1;
       end
   end
   first0 += 1;
   ωMin = ωSeq[Int(round((first0 + last0)/2))];

    # if we find an activity-disruption pair
    if (iMax != -1)&(ωMin != -1)
        brInfo1 = copy(node.brInfo);
        brInfo1[findin(pData.II,iMax)[1],findin(Ω,ωMin)[1]] = -1;
        brInfo1 = brInfoExt(pData,disData,Ω,brInfo1,ωSeq);
        brInfo2 = copy(node.brInfo);
        brInfo2[findin(pData.II,iMax)[1],findin(Ω,ωMin)[1]] = 1;
        brInfo2 = brInfoExt(pData,disData,Ω,brInfo2,ωSeq);
    end
    # add the constraints of node1 and node2
    mp1 = copy(node.mp);
    # inherit the lbCost of the branched node
    node1 = nodeType(node.lbCost,mp1,brInfo1);
    mp2 = copy(node.mp);
    # inherit the lbCost of the branched node
    node2 = nodeType(node.lbCost,mp2,brInfo2);
    return node1,node2;
end

# this is the function that identifies whether the node could not be branched anymore
function branchTest(pData,disData,Ω,node,ubCost,ϵ)
    contB = false;
    if node.lbCost <= ubCost - ϵ
        for i in pData.II
            for ω in Ω
                if node.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                    contB = true;
                end
            end
        end
    end
    return contB;
end

# this is a function that will add the unadded branching information to the master
function branchAdd(pData,disData,Ω,mp,brInfo)
    # for each non-zero item in the
    for i in pData.II
        for ω in Ω
            if brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1
                @constraint(mp,mp[:t][i] <= disData[ω].H);
            elseif brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1
                @constraint(mp,mp[:t][i] >= disData[ω].H);
            end
        end
    end
    return mp;
end

function boundCal(pData,disData,Ω,mp,lbCost,ubCost)
    mpCopy = copy(mp);
    # enforce the lower bound and the upper bound
    @constraint(mpCopy,lowerbound,pData.p0*mpCopy[:t][0] + sum(disData[ω].prDis*mpCopy[:θ][ω] for ω in Ω) >= lbCost);
    @constraint(mpCopy,upperbound,pData.p0*mpCopy[:t][0] + sum(disData[ω].prDis*mpCopy[:θ][ω] for ω in Ω) <= ubCost);

    Mub = Dict();
    Mlb = Dict();
    for i in pData.II
        @objective(mpCopy,Max,mpCopy[:t][i]);
        solve(mpCopy);
        Mub[i] = getobjectivevalue(mpCopy);
        @objective(mpCopy,Min,mpCopy[:t][i]);
        solve(mpCopy);
        Mlb[i] = getobjectivevalue(mpCopy);
    end
    return Mlb,Mub;
end

function update_brInfo(pData,disData,Ω,Mlb,Mub,brInfoT,ωSeq)
    for i in pData.II
        for ω in Ω
            if Mlb[i] >= disData[ω].H
                brInfoT[findin(pData.II,i),findin(Ω,ω)] = 1;
            end
            if Mub[i] < disData[ω].H
                brInfoT[findin(pData.II,i),findin(Ω,ω)] = -1;
            end
        end
    end
    brInfoT = brInfoExt(pData,disData,Ω,brInfoT,ωSeq);
    return brInfoT;
end
