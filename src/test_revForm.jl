using Distributed;
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,Ipopt,MathProgBase;
@everywhere using Distributions,HDF5,JLD,DelimitedFiles,Statistics,SharedArrays;
@everywhere const GUROBI_ENV = Gurobi.Env()

# test sbb
@everywhere include("header.jl");

# pathList = ["/scratch/haoxiang/current/11/",
#             "/scratch/haoxiang/current/14/",
#             "/scratch/haoxiang/current/19/",
#             "/scratch/haoxiang/current/35/",
#             "/scratch/haoxiang/current/55/",
#             "/scratch/haoxiang/current/75/"];
pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
            "/home/haoxiang/scratch/PERT_tests/current/14/",
            "/home/haoxiang/scratch/PERT_tests/current/19/",
            "/home/haoxiang/scratch/PERT_tests/current/35/",
            "/home/haoxiang/scratch/PERT_tests/current/55/",
            "/home/haoxiang/scratch/PERT_tests/current/75/"];

Ωsize = [100,200,500,1000];
sNList = [10,20,25,25];
MMList = [10,10,20,40];

fileInd = 1;
filePath = pathList[fileInd];
# compile the functions
Ωl = 1;
global Ω = 1:Ωsize[Ωl];
randNo = 1;
extBool = true;

pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1);
global pData = pData;
disDataRaw = load(pathList[fileInd]*"solData_100.jld");
disData = disDataRaw["data"][randNo];

global allSucc = findSuccAll(pData);
global distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

global sN = 20;
global MM = 5;

tFull2w,xFull2w,ubFull2w,lbFull2w,timeIter2w,treeList2w,timedecomp2w = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,true);
tFull_noMW,xFull_noMW,ubFull_noMW,lbFull_noMW,timeIter_noMW,treeList_noMW,timedecomp_noMW = partSolve_BB_para_noMW(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,21600,true);
tFull_n,xFull_n,ubFull_n,lbFull_n,timeIter_n,treeList_n,timedecomp_n = partSolve_BB_para_rev(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2);

dDict = Dict();
for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    dDict[fileInd] = Dict();
    for Ωl in 1:length(Ωsize)
        global Ω = 1:Ωsize[Ωl];
        global ϵ = 1e-2;
        global pData;
        pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize[Ωl]);
        dataRaw = load(filePath*"solData_$(Ωsize[Ωl]).jld");
        global disData = dataRaw["data"][1];

        global allSucc = findSuccAll(pData);
        global distanceDict = Dict();
        for i in pData.II
            for j in allSucc[i]
                distanceDict[i,j] = detCal(pData,i,j);
            end
        end
        # our decomposition method
        global sN = sNList[Ωl];
        global MM = MMList[Ωl];

        # A2+CS
        tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp,recordList = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,21600);
        # noMW
        tFull_noMW,xFull_noMW,ubFull_noMW,lbFull_noMW,timeIter_noMW,treeList_noMW,timedecomp_noMW = partSolve_BB_para_noMW(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,21600,true);
        # revForm
        tFull_n,xFull_n,ubFull_n,lbFull_n,timeIter_n,treeList_n,timedecomp_n = partSolve_BB_para_rev(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,21600,true);

        dDict[fileInd][Ωsize[Ωl]] = [[tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp,recordList],
                                    [tFull_noMW,xFull_noMW,ubFull_noMW,lbFull_noMW,timeIter_noMW,treeList_noMW,timedecomp_noMW],
                                    [tFull_n,xFull_n,ubFull_n,lbFull_n,timeIter_n,treeList_n,timedecomp_n]];
        save("test_Ext_revForm_1.jld","dDict",dDict);
    end
end
