addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
@everywhere const GUROBI_ENV = Gurobi.Env()

# test sbb
@everywhere include("header.jl");

pathList = ["/scratch/haoxiang/current/11/",
            "/scratch/haoxiang/current/14/",
            "/scratch/haoxiang/current/19/",
            "/scratch/haoxiang/current/35/",
            "/scratch/haoxiang/current/55/",
            "/scratch/haoxiang/current/75/"];
# pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
#             "/home/haoxiang/scratch/PERT_tests/current/14/",
#             "/home/haoxiang/scratch/PERT_tests/current/19/",
#             "/home/haoxiang/scratch/PERT_tests/current/35/",
#             "/home/haoxiang/scratch/PERT_tests/current/55/",
#             "/home/haoxiang/scratch/PERT_tests/current/75/"];

Ωsize = [100,200,500,1000,1500,2000];
sNList = [10,20,25,25,30,40];
MMList = [10,10,20,40,50,50];
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
        # Algorithm 1 without cut selection
        tFull1wo,xFull1wo,ubFull1wo,lbFull1wo,timeIter1wo,timedecomp1wo = partSolve_tightened_share(pData,disData,Ω,sN,MM,noThreads,3,5,1e-2,false);
        gapdecomp1wo = (ubFull1wo - lbFull1wo)/ubFull1wo;
        # Algorithm 1 with cut selection
        tFull1w,xFull1w,ubFull1w,lbFull1w,timeIter1w,timedecomp1w = partSolve_tightened_share(pData,disData,Ω,sN,MM,noThreads,3,5,1e-2,true);
        gapdecomp1w = (ubFull1w - lbFull1w)/ubFull1w;
        # Algorithm 2 without cut selection
        tFull2wo,xFull2wo,ubFull2wo,lbFull2wo,timeIter2wo,treeList2wo,timedecomp2wo = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,false);
        gapdecomp2wo = (ubFull2wo - lbFull2wo)/ubFull2wo;
        # Algorithm 2 with cut selection
        tFull2w,xFull2w,ubFull2w,lbFull2w,timeIter2w,treeList2w,timedecomp2w = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,10800,true);
        gapdecomp2w = (ubFull2w - lbFull2w)/ubFull2w;

        # extensive formulation
        text,xext,fext,gext,mext,timeext = extForm_cheat_new(pData,disData,Ω,sN,MM,1e-2,10800,noThreads);
        ubext = mext.objVal;
        lbext = mext.objBound;
        gapext = (ubext - lbext)/ubext;
        dDict[fileInd][Ωsize[Ωl]] = [tFull,xFull,lbFull,ubFull,gapdecomp,timedecomp,
                            text,xext,lbext,ubext,gapext,timeext];
        save("test_Ext_time.jld","dDict",dDict);
    end
end
