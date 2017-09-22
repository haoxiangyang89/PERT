# This is the main procedure of the algorithm
# for continuous crashing options and 1st case: disruption not affecting activities already started
function cutProc_Benders(pData,disData,Ω,ϵ = 1e-4)
    # initialize the first node
    nodeList = [];
    mp = createMaster(pData,disData,Ω);
    brInfo = zeros(length(pData.II),length(Ω));
    nIni = nodeType(0,mp,brInfo);
    ubCost = Inf;
    noI = 0;
    noC = 0;
    xbest = Dict();
    tbest = Dict();

    # append the first node to the nodeList
    push!(nodeList,nIni);
    ωSeq = getOmegaSeq(disData);

    # while the nodeList is not empty
    while nodeList != []
        # simply test the nodes according to the sequence they are added to the list
        # Adding other rules of node selection later !!!!!!!!!!!
        nCurrent = nodeList[1];
        # generate Benders cuts until it converges
        solve(nCurrent.mp);
        # if the current lower bound is larger than the upper bound, then prune the node
        if !(getobjectivevalue(nCurrent.mp) >= ubCost - ϵ)
            xhat = getvalue(nCurrent.mp[:x]);
            that = getvalue(nCurrent.mp[:t]);
            # update the upper bound
            ubTemp = ubCal(pData,disData,Ω,xhat,that);
            if ubTemp < ubCost
                ubCost = ubTemp;
                for i in pData.II
                    tbest[i] = that[i];
                    for j in pData.Ji[i]
                        xbest[i,j] = xhat[i,j];
                    end
                end
            end

            # update the lower bound
            if nCurrent.lbCost < getobjectivevalue(nCurrent.mp)
                lbPrev = nCurrent.lbCost;
                nCurrent.lbCost = getobjectivevalue(nCurrent.mp);
            else
                lbPrev = nCurrent.lbCost;
            end

            # stop when not improving
            while lbPrev < nCurrent.lbCost
                # plug in the xhat and that to generate a Benders cut
                for ω in Ω
                    dDω = disData[ω];
                    spCut = bGenbuild(pData,dDω,xhat,that,nCurrent.brInfo[:,findin(Ω,ω)[1]]);
                    # append the cut to the master program
                    @constraint(nCurrent.mp, nCurrent.mp[:θ][ω] >= spCut.v + sum(spCut.π[i]*(nCurrent.mp[:t][i] - that[i])
                        + sum(spCut.λ[i,j]*(nCurrent.mp[:x][i,j] - xhat[i,j]) for j in pData.Ji[i]) for i in pData.II));
                    noC += 1;
                end

                # re-solve the master problem
                solve(nCurrent.mp);
                xhat = getvalue(nCurrent.mp[:x]);
                that = getvalue(nCurrent.mp[:t]);
                # update the upper bound
                ubTemp = ubCal(pData,disData,Ω,xhat,that);
                if ubTemp < ubCost
                    ubCost = ubTemp;
                    for i in pData.II
                        tbest[i] = that[i];
                        for j in pData.Ji[i]
                            xbest[i,j] = xhat[i,j];
                        end
                    end
                end

                # update the lower bound
                if nCurrent.lbCost < getobjectivevalue(nCurrent.mp)
                    lbPrev = nCurrent.lbCost;
                    nCurrent.lbCost = getobjectivevalue(nCurrent.mp);
                else
                    lbPrev = nCurrent.lbCost;
                end
            end

            # the current node is solved to the best extent
            if branchTest(pData,disData,Ω,nCurrent,ubCost,ϵ)
                # branch the current node if it can still be branched
                # add other branching rules later !!!!!!!!!!!
                node1,node2 = branchSimple(pData,disData,Ω,nCurrent,that,ωSeq);
                if node1.brInfo != []
                    push!(nodeList,node1);
                end
                if node2.brInfo != []
                    push!(nodeList,node2);
                end
            end
        end
        # remove the current node
        shift!(nodeList);
    end
    return tbest,xbest,ubCost;
end
