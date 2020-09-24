# script to test the alternative algorithms of decomposition
using Distributed;
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,Ipopt;
@everywhere using Distributions,HDF5,JLD;
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

# altOpt: 1: no UB, no MW
#         2: UB, no MW
#         3: MW, no UB
#         4: both
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
disDataRaw = Dict();
Ωl = 6;
for fileInd in 1:4
    resultDict[fileInd] = Dict();
    disDataRaw[fileInd] = load(pathList[fileInd]*"solData_$(Ωsize[Ωl]).jld");
end
global Ω = 1:Ωsize[Ωl];
global sN = sNList[Ωl];
global MM = MMList[Ωl];

for randNo in 1:20

    for fileInd in 1:4
        filePath = pathList[fileInd];
        disData = disDataRaw[fileInd]["data"][randNo];
        pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1);
        global pData = pData;

        global allSucc = findSuccAll(pData);
        global distanceDict = Dict();
        for i in pData.II
            for j in allSucc[i]
                distanceDict[i,j] = detCal(pData,i,j);
            end
        end
        resultDict[fileInd][randNo] = [];
        for altOpt in 1:4
            if altOpt == 1
                tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_noMW(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,false);
            elseif altOpt == 2
                tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_noMW(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,true);
            elseif altOpt == 3
                tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_noUB(pData,disData,Ω,noThreads,5,6,5,1e-2,5,10800);
            elseif altOpt == 4
                tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800);
            end
            gapdecomp = (ubFull - lbFull)/ubFull;
            push!(resultDict[fileInd][randNo],[tFull,xFull,ubFull,lbFull,timedecomp]);
            save("test_alt_multi.jld","data",resultDict);
        end
    end
end

#######################################
for i in 1:4
    for j in 1:20
        if j in keys(dataalt1["data"][i])
            dataalt[i][j] = dataalt1["data"][i][j];
        elseif j in keys(dataalt2["data"][i])
            dataalt[i][j] = dataalt2["data"][i][j];
        elseif j in keys(dataalt3["data"][i])
            dataalt[i][j] = dataalt3["data"][i][j];
        elseif j in keys(dataalt4["data"][i])
            dataalt[i][j] = dataalt4["data"][i][j];
        elseif j in keys(dataalt5["data"][i])
            dataalt[i][j] = dataalt5["data"][i][j];
        end
    end
end


[dataalt[fileInd][i][1][5] for i in 1:20]
