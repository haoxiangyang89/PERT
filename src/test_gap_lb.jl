# test the gap, generate the lower bound then
using Distributed;
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
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

# read the upper bound information generated, x_hat and mean ub
bestSol = load("test_gap_ub_sol.jld");
bestSol = bestSol["sol"];

dDict = Dict();
Ωsize = [10,20,50,100,200,500,1000];
sNList = [10,20,50,100,20,25,25];
MMList = [1,1,1,1,10,20,40];
global ϵ = 1e-2;

for fileInd in 1:4
    filePath = pathList[fileInd]
    simuPath = filePath*"simuData.jld";
    data5000Raw = load(simuPath);
    dDict[fileInd] = Dict();

    for Ωl in 1:length(Ωsize)
        Ω = 1:Ωsize[Ωl];
        sN = sNList[Ωl];
        MM = MMList[Ωl];

        pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize[Ωl],20);
        dDict[fileInd][Ωsize[Ωl]] = [];
        tCanBest = bestSol[fileInd][Ωsize[Ωl]][1];
        xCanBest = bestSol[fileInd][Ωsize[Ωl]][2];
        for n in 1:20
            # randomly generate a set of disruption scenarios
            disData = disDataSet[n];
            if Ωl <= 4
                tFull,xFull,ubFull,gFull,mFull,timedecomp = extForm_cheat_new(pData,disData,Ω,sN,MM,ϵ,36000,noThreads);
                lbFull = getobjectivebound(mFull);
                gapdecomp = (ubFull - lbFull)/ubFull;
            else
                tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,ϵ,5,10800);
                gapdecomp = (ubFull - lbFull)/ubFull;
            end
            ubTemp = ubCalP(pData,disData,Ω,xCanBest,tCanBest,999999);
            push!(dDict[fileInd][Ωsize[Ωl]],[tFull,xFull,lbFull,ubFull,gapdecomp,ubTemp]);
            save("test_gap.jld","lbDict",dDict);
        end
    end
end
