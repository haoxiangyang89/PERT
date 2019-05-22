addprocs(30);
global noThreads = 30;
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
# pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
#             "/home/haoxiang/scratch/PERT_tests/current/14/",
#             "/home/haoxiang/scratch/PERT_tests/current/19/",
#             "/home/haoxiang/scratch/PERT_tests/current/35/",
#             "/home/haoxiang/scratch/PERT_tests/current/55/",
#             "/home/haoxiang/scratch/PERT_tests/current/75/"];

Ωsize = [100,200,500,1000,1500,2000];
sNList = [10,20,20,20,30,40];
MMList = [10,10,25,50,50,50];
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
        tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,1,1e-2,5,1000);
        gapdecomp = (ubFull - lbFull)/ubFull;

        # extensive formulation
        tic();
        text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,1e-4,10800);
        timeext = toc();
        ubmp = mext.objVal;
        lbmp = mext.objBound;
        gapext = (mext.objVal - mext.objBound)/mext.objVal;
        ubext = ubCal(pData,disData,Ω,xext,text,9999999);
        dDict[fileInd][Ωsize[Ωl]] = [tFull,xFull,lbFull,ubFull,gapdecomp,timedecomp,
                            text,xext,lbmp,ubmp,ubext,gapext,timeext];
        save("test_Ext_time.jld","dDict",dDict);
    end
end
