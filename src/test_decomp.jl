# test the decomposition for a single problem
addprocs(31);
global noThreads = 31;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

pathList = ["/scratch/haoxiang/current/11/",
            "/scratch/haoxiang/current/14/",
            "/scratch/haoxiang/current/19/",
            "/scratch/haoxiang/current/35/",
            "/scratch/haoxiang/current/55/",
            "/scratch/haoxiang/current/75/"];

dDict = Dict();
ubDict = Dict();
Ωsize = 500;
global Ω = 1:Ωsize;
global ϵ = 1e-2;
global sN = 2;
global MM = 25;

fileInd = 5;
filePath = pathList[fileInd];
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize);
global pData = pData;
dataRaw = load(filePath*"solData_500.jld");
disData = dataRaw["data"][1];
# global disData = disDataSet[1];
dDict[fileInd] = [];
ubDict[fileInd] = [];

global allSucc = findSuccAll(pData);
global distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

tic();
#include("partSolve_Callback_tightened_sol.jl");
tFull,xFull,ubFull,lbFull,timeIter,treeList = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,1,1e-2,5,1000,21600);
timeFull = toc();

save("test_decomp_55.jld","data",[tFull,xFull,ubFull,lbFull,timeIter,treeList,timeFull]);
