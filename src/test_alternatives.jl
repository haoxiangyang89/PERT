# script to test the alternative algorithms of decomposition
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
@everywhere const GUROBI_ENV = Gurobi.Env()

# test sbb
@everywhere include("header.jl");

# pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
#             "/home/haoxiang/scratch/PERT_tests/current/14/",
#             "/home/haoxiang/scratch/PERT_tests/current/19/",
#             "/home/haoxiang/scratch/PERT_tests/current/35/",
#             "/home/haoxiang/scratch/PERT_tests/current/55/",
#             "/home/haoxiang/scratch/PERT_tests/current/75/"];
pathList = ["/scratch/haoxiang/current/11/",
               "/scratch/haoxiang/current/14/",
               "/scratch/haoxiang/current/19/",
               "/scratch/haoxiang/current/35/",
               "/scratch/haoxiang/current/55/",
               "/scratch/haoxiang/current/75/"];

# altOpt: 1: no UB, no MW
#         2: UB, no MW
#         3: MW, no UB
#         4: both
altOpt = 1
Ωsize = [10,20,50,100,200,500,1000,2000];
sNList = [10,0,0,0,20,25,20,40];
MMList = [1,0,0,0,10,20,50,50];

fileInd = 1;
filePath = pathList[fileInd];
# compile the functions
Ωl = 1;
global Ω = 1:Ωsize[Ωl];
randNo = 1;
extBool = true;

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

tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_noMW(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,true);
tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_noUB(pData,disData,Ω,noThreads,5,6,5,1e-2,5,10800);
tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800);


resultDict = Dict();

for fileInd in 1:4
    Ωl = 6;
    filePath = pathList[fileInd];
    global Ω = 1:Ωsize[Ωl];
    randNo = 1;

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
    resultDict[fileInd] = [];

    for altOpt in 1:4
        if altOpt == 1
            tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_noMW(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,false);
        elseif altOpt == 2
            tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_noMW(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,true);
        elseif altOpt == 3
            tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_noUB(pData,disData,Ω,noThreads,5,6,5,1e-2,5,10800);
        elseif altOpt == 4
            tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,false);
        end
        gapdecomp = (ubFull - lbFull)/ubFull;
        push!(resultDict[fileInd],[tFull,xFull,ubFull,lbFull,timedecomp]);
        save("test_alt.jld","data",resultDict);
    end
end
