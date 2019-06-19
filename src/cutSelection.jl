# cut selection for the decomposition method
function examineCuts_count(disData,Ω,cutSel,cutSet,that,xhat,θhat,yhat,cutThreshold)
    # check for each cut whether it is tight, if not update the counts
    cutyn = [];
    for nc in 1:length(cutSet)
        # how many rounds have been through
        for l in 1:length(cutSet[nc][2])
            # for each cut
            ω = cutSet[nc][2][l][1]
            cutV = cutSet[nc][2][l][2];
            for i in pData.II
                cutV += cutSet[nc][2][l][3][i]*(that[i] - cutSet[nc][1][1][i]);
                for j in pData.Ji[i]
                    cutV += cutSet[nc][2][l][4][i,j]*(xhat[i,j] - cutSet[nc][1][2][i,j]);
                end
                for par in 1:length(cutSet[nc][1][4][i])
                    cutV += cutSet[nc][2][l][5][i,par]*(sum(yhat[i,parNew] for parNew in 1:length(divSet[i]) if revPar(cutSet[nc][1][4][i],divSet[i][parNew]) == par) - cutSet[nc][1][3][i,par]);
                end
            end
            if abs(θhat[ω] - cutV)/θhat[ω] > 1e-4
                # not tight
                cutSel[nc,l] += 1;
            else
                cutSel[nc,l] = 0;
            end
        end
    end
    return cutSel;
end

function selectCuts(cutSel,cutSet,cutThreshold)
    cutSelNew = Dict();
    for nc in 1:length(cutSet)
        leftSet = [l for l in 1:length(cutSet[nc][2]) if cutSel[nc,l] <= cutThreshold];
        cutSet[nc][2] = cutSet[nc][2][leftSet];
        for l1 in 1:length(leftSet)
            cutSelNew[nc,l1] = cutSel[nc,leftSet[l1]];
        end
    end
    return cutSelNew,cutSet;
end

function examineCuts_count_2(pData,disData,Ω,cutSet,divSet,that,xhat,θhat,yhat)
    # check for each cut whether it is tight, if not update the counts
    cutSel = [];
    for nc in 1:length(cutSet)
        # how many rounds have been through
        for l in 1:length(cutSet[nc][2])
            # for each cut
            ω = cutSet[nc][2][l][1]
            cutV = cutSet[nc][2][l][2];
            for i in pData.II
                cutV += cutSet[nc][2][l][3][i]*that[i];
                for j in pData.Ji[i]
                    cutV += cutSet[nc][2][l][4][i,j]*xhat[i,j];
                end
                for par in 1:length(cutSet[nc][1][i])
                    cutV += cutSet[nc][2][l][5][i,par]*(sum(yhat[i,parNew] for parNew in 1:length(divSet[i]) if revPar(cutSet[nc][1][i],divSet[i][parNew]) == par));
                end
            end
            if abs(θhat[ω] - cutV)/θhat[ω] <= 1e-4
                # tight
                push!(cutSel,(nc,l));
            end
        end
    end
    return cutSel;
end

function examineCuts_count_3(pData,disData,Ω,cutSet,divSet,that,xhat,θhat,yhat,IJPair,IPPair)
    # check for each cut whether it is tight, if not update the counts
    cutSel = Dict();
    for nc in 1:length(cutSet)
        cutSel[nc] = [];
        # how many rounds have been through
        for l in 1:length(cutSet[nc][1])
            # for each cut
            ω = cutSet[nc][1][l];
            cutV = cutSet[nc][5][l];
            for i in pData.II
                cutV += cutSet[nc][2][findfirst(pData.II,i),l]*that[i];
                for j in pData.Ji[i]
                    cutV += cutSet[nc][3][findfirst(IJPair,(i,j)),l]*xhat[i,j];
                end
                for par in 1:length(divSet[i])
                    cutV += cutSet[nc][4][findfirst(IPPair,(i,par)),l]*yhat[i,par];
                end
            end
            if abs(θhat[ω] - cutV)/θhat[ω] <= 1e-4
                # tight
                push!(cutSel[nc],l);
            end
        end
    end
    return cutSel;
end

function selectCuts2(cutSet,cutSetSel)
    cutSetNew = [];
    for nc in 1:length(cutSet)
        leftSet = [l for l in 1:length(cutSet[nc][2]) if (nc,l) in cutSetSel];
        push!(cutSetNew, [cutSet[nc][1],cutSet[nc][2][leftSet]]);
    end
    return cutSetNew;
end

function selectCuts3(cutSet,cutSel)
    cutSetNew = deepcopy(cutSet);
    for nc in 1:length(cutSet)
        if cutSet[nc][1] != []
            # if the original set is not empty
            cutSetNew[nc][1] = cutSet[nc][1][cutSel[nc]];
            cutSetNew[nc][2] = cutSet[nc][2][:,cutSel[nc]];
            cutSetNew[nc][3] = cutSet[nc][3][:,cutSel[nc]];
            cutSetNew[nc][4] = cutSet[nc][4][:,cutSel[nc]];
            cutSetNew[nc][5] = cutSet[nc][5][cutSel[nc]];
        end
    end
    return cutSetNew;
end
