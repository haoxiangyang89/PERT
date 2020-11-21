# compare the output of A2 and Gurobi extensive formulation
using Distributed;
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,Ipopt,MathProgBase;
@everywhere using Distributions,HDF5,JLD,DelimitedFiles,Statistics,SharedArrays;
@everywhere const GUROBI_ENV = Gurobi.Env();

# test sbb
@everywhere include("header.jl");

pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
            "/home/haoxiang/scratch/PERT_tests/current/14/",
            "/home/haoxiang/scratch/PERT_tests/current/19/",
            "/home/haoxiang/scratch/PERT_tests/current/35/",
            "/home/haoxiang/scratch/PERT_tests/current/55/",
            "/home/haoxiang/scratch/PERT_tests/current/75/"];
# pathList = ["/scratch/haoxiang/current/11/",
#                "/scratch/haoxiang/current/14/",
#                "/scratch/haoxiang/current/19/",
#                "/scratch/haoxiang/current/35/",
#                "/scratch/haoxiang/current/55/",
#                "/scratch/haoxiang/current/75/"];

altOpt = 1
Ωsize = [10,20,50,100,200,500,1000,2000];
sNList = [10,0,0,0,20,25,25,40];
MMList = [1,0,0,0,10,20,40,50];

fileInd = 1;
filePath = pathList[fileInd];
# compile the functions
Ωl = 5;
global Ω = 1:Ωsize[Ωl];
randNo = 1;
extBool = true;
global ϵ = 1e-2;

pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1);
global pData = pData;
disDataRaw = load(pathList[fileInd]*"solData_$(Ωsize[Ωl]).jld");
disData = disDataRaw["data"][randNo];

global allSucc = findSuccAll(pData);
global distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

global sN = sNList[Ωl];
global MM = MMList[Ωl];

# prerun to compile the codes
tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp,recordList = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800);
text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω[1:5],ϵ,999999);

#######################################################################################################
# run the full-scale experiments
fileInd = 2;
Ωl = 7;
global Ω = 1:Ωsize[Ωl];

pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1);
global pData = pData;
disDataRaw = load(pathList[fileInd]*"solData_$(Ωsize[Ωl]).jld");
disData = disDataRaw["data"][randNo];
global allSucc = findSuccAll(pData);
global distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

global sN = sNList[Ωl];
global MM = MMList[Ωl];
tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp,recordList = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,21600);
text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,ϵ,999999,30,"progress_Gurobi.log");
save("test_progressComp.jld","data",[tFull,xFull,ubFull,lbFull,timedecomp,recordList]);
