function cutValue(pData,disData,text,xext,gext,θext,cutSet)
    H = Dict();
    H[0] = 0;
    H[length(Ω)+1] = Tmax;
    for ω in Ω
        H[ω] = disData[ω].H;
    end

    θs = zeros(length(Ω),length(cutSet));
    # for each batch of cut
    for nc in 1:length(cutSet)
        # identify the value of y variables
        divSet = cutSet[nc][1][4];
        yext = Dict();
        for i in pData.II
            for par in 1:length(divSet[i])
                if (text[i] > H[divSet[i][par].startH])&(text[i] < H[divSet[i][par].endH])
                    yext[i,par] = 1;
                elseif text[i] == H[divSet[i][par].startH]
                    if divSet[i][par].startH > 0
                        if gext[i,divSet[i][par].startH] == 0
                            yext[i,par] = 0;
                        else
                            yext[i,par] = 1;
                        end
                    else
                        yext[i,par] = 1;
                    end
                elseif text[i] == H[divSet[i][par].endH]
                    if divSet[i][par].startH < length(Ω) + 1
                        if gext[i,divSet[i][par].endH] == 1
                            yext[i,par] = 0;
                        else
                            yext[i,par] = 1;
                        end
                    else
                        yext[i,par] = 1;
                    end
                else
                    yext[i,par] = 0;
                end
            end
        end
        # for each scenario
        for ω in Ω
            if cutSet[nc][2][ω] != []
                θs[ω,nc] = cutSet[nc][2][ω][1] + sum(cutSet[nc][2][ω][2][i]*(text[i] - cutSet[nc][1][1][i]) for i in pData.II) +
                    sum(sum(cutSet[nc][2][ω][3][i,j]*(xext[i,j] - cutSet[nc][1][2][i,j]) for j in pData.Ji[i]) for i in pData.II if i != 0) +
                    sum(sum(cutSet[nc][2][ω][4][i,par]*(yext[i,par] - cutSet[nc][1][3][i,par]) for par in 1:length(cutSet[nc][1][4][i])) for i in pData.II);
            else
                θs[ω,nc] = 0;
            end
        end
    end

    for nc in 1:length(cutSet)
        for ω in Ω
            if θs[ω,nc] > θext[ω]
                print((ω,nc)," ",θs[ω,nc]," ",θext[ω]);
            end
        end
    end
    return θs;
end

function cutValueY(pData,disData,text,xext,yext,θext,cutSet,divSet)
    H = Dict();
    H[0] = 0;
    H[length(Ω)+1] = Tmax;
    for ω in Ω
        H[ω] = disData[ω].H;
    end

    θs = Dict();
    # for each batch of cut
    for nc in 1:length(cutSet)
        # for each scenario
        for l in 1:length(cutSet[nc][2])
            ω = cutSet[nc][2][l][1];
            vk = cutSet[nc][2][l][2];
            πk = cutSet[nc][2][l][3];
            λk = cutSet[nc][2][l][4];
            γk = cutSet[nc][2][l][5];
            θs[ω,nc] = vk + sum(πk[i]*(text[i] - cutSet[nc][1][1][i]) for i in pData.II);
            for i in pData.II
                for j in pData.Ji[i]
                    θs[ω,nc] += λk[i,j]*(xext[i,j] - cutSet[nc][1][2][i,j]);
                end
                for par in 1:length(cutSet[nc][1][4][i])
                    ysum = 0;
                    for parNew in 1:length(divSet[i])
                        if revPar(cutSet[nc][1][4][i],divSet[i][parNew]) == par
                            ysum += yext[i,parNew];
                        end
                    end
                    θs[ω,nc] += γk[i,par]*(ysum - cutSet[nc][1][3][i,par]);
                end
            end
        end
    end

    for (ω,nc) in keys(θs)
        if θs[ω,nc] > θext[ω]
            println((ω,nc)," ",θs[ω,nc]," ",θext[ω]);
        end
    end
    return θs;
end
