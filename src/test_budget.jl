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


# filePath = "/Users/haoxiangyang/Desktop/PERT_tests/14_Lognormal_Exponential/"
dDict = Dict();
Ωsize = [10,20,50,100,200,300,400,500,750,1000,1500,2000];
sNList = [0,0,0,0,20,20,20,20,25,20,30,40];
MMList = [0,0,0,0,10,15,20,25,30,50,50,50];
for fileInd in 1:length(pathList)
    filePath = pathList[fileInd]
    simuPath = filePath*"simuData.jld";
    data5000Raw = load(simuPath);

    dDict[fileInd] = Dict();
    disData1 = data5000Raw["data"][1];
    Ω1 = 1:length(disData1);
    pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1);
    global pData = pData;
    global allSucc = findSuccAll(pData);
    global distanceDict = Dict();
    for i in pData.II
        for j in allSucc[i]
            distanceDict[i,j] = detCal(pData,i,j);
        end
    end

    for Ωl in 1:length(Ωsize)
        global Ω = 1:Ωsize[Ωl];
        disDataRaw = load(pathList[fileInd]*"solData_$(Ωsize[Ωl]).jld");
        disDataSet = disDataRaw["data"];
        global ϵ = 1e-2;
        dDict[fileInd][Ωsize[Ωl]] = [];
        ubbudget = [];
        n = 1;
        while n <= 20
            # try
            global disData = disDataSet[n];

            # our decomposition method
            if Ωl <= 4
                tic();
                tFull,xFull,ubFull,gFull,mFull = extForm_cheat(pData,disData,Ω,1e-4,999999,noThreads);
                timedecomp = toc();
                gapdecomp = 0;
                lbFull = getobjectivebound(mFull);
            else
                global sN = sNList[Ωl];
                global MM = MMList[Ωl];
                global nSplit = 5;
                tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,1,1e-2,5,1000);
                gapdecomp = (ubFull - lbFull)/ubFull;
            end

            ubTemp = ubCalP(pData,disData1,Ω1,xFull,tFull,999999);
            push!(dDict[fileInd][Ωsize[Ωl]],[tFull,xFull,lbFull,ubFull,gapdecomp,timedecomp,ubTemp]);
            save("test_Ext_budget.jld","dDict",dDict);
            println("========= Case $(fileInd) Budget $(Ωsize[Ωl]) Sample $(n) Processed, Lower bound $(lbFull), Upper bound $(ubFull) =========");
            n += 1;
            # catch
            # println("Error in Data!");
            # push!(ErrorData,(fileInd,disData));
            # end
        end
    end
end
