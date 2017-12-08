# pull binary variables to front and manage a branch and bound tree

function branchPull(pData,disData,Ω,node,lb,ub)
    ωSeq,ωDict = getOmegaSeq(disData);
    # obtain the number of zeros in brInfo for each activity
    # select which activity is branched on
    noZeros = Dict();
    iMax = -1;
    maxZero = 0;
    for i in pData.II
        noZeros[i] = sum([nCurrent.brInfo[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0 for ω in Ω]);
        if noZeros[i] > maxZero
            maxZero = noZeros[i];
            iMax = i;
        end
    end
    # iMax = -1;
    # maxRange = 0;
    # for i in pData.II
    #     println(i," ",node.tmaxD[i] - node.tminD[i]);
    #     if node.tmaxD[i] - node.tminD[i] > maxRange
    #         maxRange = node.tmaxD[i] - node.tminD[i];
    #         iMax = i;
    #     end
    # end

    # select which scenario's time we branch against
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
    # ωMin = -1;
    # minRange = 999999;
    # for ω in Ω
    #     if abs(that[iMax] - disData[ω].H) < minRange
    #         minRange = abs(that[iMax] - disData[ω].H);
    #         ωMin = ω;
    #     end
    # end

    # branch!
    if (iMax != -1)&(ωMin != -1)
        # branch down
        tMaxD1 = copy(node.tmaxD);
        tMinD1 = copy(node.tminD);
        brInfo1 = copy(node.brInfo);
        tMaxD1[iMax] = disData[ωMin].H - 1e-10;
        bp = buildTighten(pData,disData,Ω,node.cutSet,ub,lb,tMaxD1,tMinD1,node.brInfo);
        maxt = Dict();
        mint = Dict();
        node1state = true;
        for i in pData.II
            # change the objective function value
            @objective(bp, Max, bp[:t][i]);
            bpStatus = solve(bp);
            if bpStatus == :Optimal
                maxt[i] = getobjectivevalue(bp);
                tMaxD1[i] = min(maxt[i],tMaxD1[i]);
                @objective(bp, Min, bp[:t][i]);
                solve(bp);
                mint[i] = getobjectivevalue(bp);
                tMinD1[i] = max(mint[i],tMinD1[i]);
                for ω in Ω
                    if maxt[i] < disData[ω].H
                        if brInfo1[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                            brInfo1[findin(pData.II,i)[1],findin(Ω,ω)[1]] = -1;
                        end
                    end
                    if mint[i] >= disData[ω].H
                        if brInfo1[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                            brInfo1[findin(pData.II,i)[1],findin(Ω,ω)[1]] = 1;
                        end
                    end
                end
            else
                node1state = false;
            end
        end
        brInfo1 = brInfoExt(pData,disData,Ω,brInfo1,ωSeq);
        node1 = nodeTypeP(lb,brInfo1,tMaxD1,tMinD1,node.cutSet,node1state);

        # branch up
        tMaxD2 = copy(node.tmaxD);
        tMinD2 = copy(node.tminD);
        brInfo2 = copy(node.brInfo);
        tMinD2[iMax] = disData[ωMin].H;
        bp = buildTighten(pData,disData,Ω,node.cutSet,ub,lb,tMaxD2,tMinD2,node.brInfo);
        maxt = Dict();
        mint = Dict();
        node2state = true;
        for i in pData.II
            # change the objective function value
            @objective(bp, Max, bp[:t][i]);
            bpStatus = solve(bp);
            if bpStatus == :Optimal
                maxt[i] = getobjectivevalue(bp);
                tMaxD2[i] = min(maxt[i],tMaxD2[i]);
                @objective(bp, Min, bp[:t][i]);
                solve(bp);
                mint[i] = getobjectivevalue(bp);
                tMinD2[i] = max(mint[i],tMinD2[i]);
                for ω in Ω
                    if maxt[i] < disData[ω].H
                        if brInfo2[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                            brInfo2[findin(pData.II,i)[1],findin(Ω,ω)[1]] = -1;
                        end
                    end
                    if mint[i] >= disData[ω].H
                        if brInfo2[findin(pData.II,i)[1],findin(Ω,ω)[1]] == 0
                            brInfo2[findin(pData.II,i)[1],findin(Ω,ω)[1]] = 1;
                        end
                    end
                end
            else
                node2state = false;
            end
        end
        brInfo2 = brInfoExt(pData,disData,Ω,brInfo2,ωSeq);
        node2 = nodeTypeP(lb,brInfo2,tMaxD2,tMinD2,node.cutSet,node2state);
    end
    # return the branched nodes
    return node1,node2;
end
