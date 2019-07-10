# illustration example for the branch-and-bound algorithm
addprocs(3);
global noThreads = 3;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
@everywhere const GUROBI_ENV = Gurobi.Env();
# test sbb
@everywhere include("header.jl");

filePath = "/Users/haoxiangyang/Desktop/PERT_tests/current/5/"
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,4);
global pData = pData;
global allSucc = findSuccAll(pData);
global distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

Ω = 1:4;
disData = disDataSet[1];
for ω in Ω
    disData[ω].H = ω;
end
global ϵ = 1e-2;
tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_noMW(pData,disData,Ω,noThreads,1,1,1e-4,1,2000);
