# cut selection for the decomposition method
function examineCuts_count(disData,Ω,cutSel,cutSet,that,xhat,θhat,yhat,cutThreshold)
    # check for each cut whether it is tight, if not update the counts
    cutyn = [];
    for nc in 1:length(cutSet)
        # how many rounds have been through
        for ω in Ω
            # for each scenario
            if cutSet[nc][2][ω] != []
                # there is a cut
                cutV = cutSet[nc][2][ω][1];
                for i in pData.II
                    cutV += cutSet[nc][2][ω][2][i]*(that[i] - cutSet[nc][1][1][i]);
                    for j in pData.Ji[i]
                        cutV += cutSet[nc][2][ω][3][i,j]*(xhat[i,j] - cutSet[nc][1][2][i,j]);
                    end
                    for par in 1:length(cutSet[nc][1][4][i])
                        cutV += cutSet[nc][2][ω][4][i,par]*(sum(yhat[i,parNew] for parNew in 1:length(divSet[i]) if revPar(cutSet[nc][1][4][i],divSet[i][parNew]) == par) - cutSet[nc][1][3][i,par]);
                    end
                end
                if θhat[ω] > cutV + 1e-4
                    # not tight
                    cutSel[nc,ω] += 1;
                else
                    cutSel[nc,ω] = 0;
                end
            end
            if (cutSel[nc,ω] > cutThreshold)&(!((nc,ω) in cutyn))
                push!(cutyn,(nc,ω));
            end
        end
    end
    return cutSel,cutyn;
end

function examineCuts_slack(disData,Ω,cutSlack,cutSet,that,xhat,θhat,yhat)
    # check for each cut whether it is tight, if not update the counts

    return cutSlack;
end
