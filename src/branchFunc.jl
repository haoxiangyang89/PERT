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
        # add the branching constraint at each iteration
        # for i in pData.II
        #     for ω in Ω
        #         if brInfo1[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1
        #             @constraint(mp1,mp1[:t][i] <= disData[ω].H - 1e-6);
        #         elseif brInfo1[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1
        #             @constraint(mp1,mp1[:t][i] >= disData[ω].H);
        #         end
        #     end
        # end
        # inherit the lbCost of the branched node
        node1 = nodeType(node.lbCost,mp1,brInfo1);
    end
    mp2 = copy(node.mp);
    if brInfo2 != []
        # for i in pData.II
        #     for ω in Ω
        #         if brInfo2[findin(pData.II,i)[1],findin(Ω,ω)[1]] == -1
        #             @constraint(mp2,mp2[:t][i] <= disData[ω].H - 1e-6);
        #         elseif brInfo2[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1
        #             @constraint(mp2,mp2[:t][i] >= disData[ω].H);
        #         end
        #     end
        # end
        # inherit the lbCost of the branched node
        node2 = nodeType(node.lbCost,mp2,brInfo2);
    end
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
                @constraint(mp,mp[:t][i] <= disData[ω].H - 1e-6);
            elseif brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 1
                @constraint(mp,mp[:t][i] >= disData[ω].H);
            end
        end
    end
    return mp;
end
