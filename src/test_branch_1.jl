addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

pathList = ["/home/haoxiang/PERT_tests/current/11/",
            "/home/haoxiang/PERT_tests/current/14/",
            "/home/haoxiang/PERT_tests/current/19/",
            "/home/haoxiang/PERT_tests/current/35/",
            "/home/haoxiang/PERT_tests/current/55/",
            "/home/haoxiang/PERT_tests/current/75/"];

# pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
#             "/home/haoxiang/scratch/PERT_tests/current/14/",
#             "/home/haoxiang/scratch/PERT_tests/current/19/",
#             "/home/haoxiang/scratch/PERT_tests/current/35/",
#             "/home/haoxiang/scratch/PERT_tests/current/55/",
#             "/home/haoxiang/scratch/PERT_tests/current/75/"];

fileInd = 2;
filePath = pathList[fileInd];
Ωsize = 500;
global Ω = 1:Ωsize;
global ϵ = 1e-2;
global pData;
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize);
global pData = pData;
dataRaw = load("14_test1.jld");
global disData = dataRaw["disData"];

global allSucc = findSuccAll(pData);
global distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

global sN = 25;
global MM = 20;
global nSplit = 5;
data = [];

global bAlt = 1;
tempTimer = time();
include("partSolve_Callback_tightened_sol.jl");
timeFull = time() - tempTimer;
gapFull = (ubCost - lbCost)/ubCost;
ubFull = ubCost;
lbFull = lbCost;
xFull = deepcopy(xbest);
tFull = deepcopy(tbest);
push!(data,[ubFull,tFull,xFull,lbFull,timeFull]);
save("14_test1_Results_1.jld","data",data);
