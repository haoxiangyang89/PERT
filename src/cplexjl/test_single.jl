# test the running time given a number of scenarios and obtain the result

addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
            "/home/haoxiang/scratch/PERT_tests/current/14/",
            "/home/haoxiang/scratch/PERT_tests/current/19/",
            "/home/haoxiang/scratch/PERT_tests/current/35/",
            "/home/haoxiang/scratch/PERT_tests/current/55/",
            "/home/haoxiang/scratch/PERT_tests/current/75/"];

# pathList = ["/Users/haoxiangyang/Desktop/PERT_tests/current/11/",
#             "/Users/haoxiangyang/Desktop/PERT_tests/current/14/",
#             "/Users/haoxiangyang/Desktop/PERT_tests/current/19/",
#             "/Users/haoxiangyang/Desktop/PERT_tests/current/35/",
#             "/Users/haoxiangyang/Desktop/PERT_tests/current/55/",
#             "/Users/haoxiangyang/Desktop/PERT_tests/current/75/"];

#fileInd = 1;
Ωsize = [10,20,50,100,200,500,1000,2000];
sNList = [0,0,0,0,20,25,20,40];
MMList = [0,0,0,0,10,20,50,50];
for fileInd in 1:4
    filePath = pathList[fileInd];

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
    if extBool
        tic();
        #tFull,xFull,ubFull,gFull,mFull = extForm_cheat(pData,disData,Ω,1e-4,999999,noThreads);
        tFull,xFull,ubFull,gFull,mFull = extForm_cheat_new(pData,disData,Ω,sN,MM,1e-2,999999,noThreads);
        timedecomp = toc();
        lbFull = getobjectivebound(mFull);
        gapdecomp = (ubFull - lbFull)/ubFull;
    else
        tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,1e-2,5,10800,5);
        gapdecomp = (ubFull - lbFull)/ubFull;
    end

    save("test_single_$(fileInd)_$(Ωsize[Ωl])_$(extBool).jld","data",[tFull,xFull,ubFull,lbFull,timedecomp,gapdecomp]);
end
