function addCuts(pData,Ω,mp,πdict,γdict,λdict,vk,that,xhat,Ghat,cutSet,partRev)
    # add the linear Benders cuts
    for ω in Ω
        Ghatω = Ghat[ω];
        @constraint(mp,mp[:θ][ω] >= vk[ω] + sum(πdict[ω][i]*(mp[:t][i] - that[i]) for i in pData.II) +
            sum(sum(λdict[ω][i,j]*(mp[:x][i,j] - xhat[i,j]) for j in pData.Ji[i]) for i in pData.II) +
            sum(γdict[ω][i]*(mp[:G][i,partRev[i][ω]] - Ghatω[i]) for i in pData.II));
        push!(cutSet[ω],(πdict[ω],γdict[ω],λdict[ω],vk[ω],that,xhat,Ghatω));
    end
    return mp,cutSet;
end
