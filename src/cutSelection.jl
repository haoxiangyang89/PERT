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

function examineCuts_slack(disData,Ω,cutSlack,cutSet,that,xhat,θhat,yhat)
    # check for each cut whether it is tight, if not update the counts

    return cutSlack;
end
