##########################################################################################
# obtain the RHS of the cut
cutSet[ω][nc][4] + sum(cutSet[ω][nc][1][i]*(that[i] - cutSet[ω][nc][5][i]) for i in pData.II) +
            sum(sum(cutSet[ω][nc][2][i,j]*(xhat[i,j] - cutSet[ω][nc][6][i,j]) for j in pData.Ji[i]) for i in pData.II if i != 0) +
            sum(cutSet[ω][nc][3][i,ω]*(Ghat[i,ω] - cutSet[ω][nc][7][i,ω]) for i in pData.II if (i,ω) in keys(cutSet[ω][nc][7]));

##########################################################################################
# obtain the 50/500 test data
pInputAdd = "test_14_P.csv";
kInputAdd = "test_14_K.csv";
ΩInputAdd = "test_14_Omega_full.csv";
ϕInputAdd = "test_14_Phi_full.csv";
pData = readInP(pInputAdd,kInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
# disData,Ω = autoUGen("LogNormal",[log(35),0.5],nameD,dparams,50,1 - pData.p0);
disData,Ω = autoUGen("LogNormal",[log(35),0.5],nameD,dparams,500,1 - pData.p0);
disData = orderdisData(disData,Ω);
text,xext,fext,gext,mp = extForm_cheat(pData,disData,Ω);
tBest = Dict();
xBest = Dict();
GBest = Dict();
for i in pData.II
    tBest[i] = text[i];
    for j in pData.Ji[i]
        xBest[i,j] = xext[i,j];
    end
    for ω in Ω
        GBest[i,ω] = gext[i,ω];
    end
end
save("testData_500.jld","text",tBest,"xext",xBest,"gext",GBest,"fext",fext,"disData",disData);

##########################################################################################
θext = Dict();
for ω in Ω
    θext[ω] = subIntG(pData,disData[ω],xext,text,GextDict[ω]);
end

##########################################################################################
# Within the iteration, check the cut validity
θInt = Dict();
for ω in Ω
    θInt[ω] = subIntGMixed(pData,disData[ω],xhat,that,ω,Ghat);
    println(vk[ω]," ",θInt[ω]);
end
for ω in Ω
    if vk[ω] > θInt[ω] + 1e-6
        println(ω);
    end
end
for ω in Ω
    if vk[ω] > θext[ω] + 1e-6
        println(ω);
    end
end
for ω in Ω
    if θlb[ω] > θext[ω] + 1e-6
        println(ω);
    end
end
for ω in Ω
    if subMixedDict[ω] > θext[ω] + 1e-6
        println(ω," ",subMixedDict[ω]," ",θext[ω]," ");
    end
end
