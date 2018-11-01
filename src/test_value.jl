addprocs(20);
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

filePath = "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/";
Ωsize = 500;
Ω = 1:Ωsize;
ϵ = 1e-2;
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize);
disData = disDataSet[1];

allSucc = findSuccAll(pData);
distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end
#disData = disDataSet[1];

# deterministic solution
tdet,xdet,fdet = detBuild(pData);
ubdet = ubCalP(pData,disData,Ω,xdet,tdet,999999);

# expected solution
eH = mean(buildDistrn(nameH,Hparams));
ed = Dict();
for i in pData.II
    if i == 0
        ed[i] = 0;
    else
        ed[i] = mean(buildDistrn(nameD,dparams[i]));
    end
end
texp,xexp,fexp,Gexp,mexp = expModel(pData,eH,ed);
ubexp = ubCalP(pData,disData,Ω,xexp,texp,999999);

# our decomposition method
tic();
#tbest,xbest,lbCost,ubCost = partitionSolve(pData,disData,0.001);
#include("partSolve_Callback.jl");
include("partSolve_Callback_tightened.jl");
timedecomp = toc();
gapdecomp = (ubCost - lbCost)/ubCost;
