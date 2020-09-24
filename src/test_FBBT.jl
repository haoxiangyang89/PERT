# test the effectiveness of FBBT
using Distributed;
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

Ωsize = [10,20,50,100,200,500,1000,2000];
sNList = [10,0,0,0,20,25,25,40];
MMList = [1,0,0,0,10,20,40,50];

# precompile
fileInd = 1;
filePath = pathList[fileInd];
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
tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800);

data = Dict();
for fileInd in 1:4
    Ωl = 7;
    filePath = pathList[fileInd];
    data[fileInd] = Dict();
    global Ω = 1:Ωsize[Ωl];
    for randNo in 1:20
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

        # test the solution process with/out FBBT, initial LB is recorded in the treeList
        # with
        tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,true);
        # without
        tFullo,xFullo,ubFullo,lbFullo,timeItero,treeListo,timedecompo = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,true,false);
        data[fileInd][randNo] = [tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp,
                        tFullo,xFullo,ubFullo,lbFullo,timeItero,treeListo,timedecompo];
        save("test_FBBT.jld","data",data);
    end
end
