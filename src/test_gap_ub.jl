# test the gap, generate the upper bound first
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

ubDict = Dict();
Ωsize = [10,20,50,100,200,500,1000];
sNList = [10,20,50,100,20,25,25];
MMList = [1,1,1,1,10,20,40];
data5000RawD = load("data5000Raw_tot.jld");
data5000RawD = data5000RawD["data"];
global ϵ = 1e-2;
for fileInd in 1:4
    ubDict[fileInd] = Dict();
end

for Ωl in 1:length(Ωsize)
    disDataRaw = Dict();
    for fileInd in 1:4
        disDataRaw[fileInd] = load(pathList[fileInd]*"solData_$(Ωsize[Ωl]).jld");
        ubDict[fileInd][Ωsize[Ωl]] = Dict();
    end
    Ω = 1:Ωsize[Ωl];
    sN = sNList[Ωl];
    MM = MMList[Ωl];

    for randNo in 1:20

        for fileInd in 1:4
            filePath = pathList[fileInd]
            disData1 = data5000RawD[fileInd];
            Ω1 = 1:length(disData1);

            # generate a candidate solution by SAA
            pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize[Ωl]);
            disData = disDataSet[1];
            tCan,xCan,ubCan,lbCan,timeIterCan,treeListCan,timedecompCan = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,ϵ,5,10800);
            ub5000 = ubCalP(pData,disData1,Ω1,xCan,tCan,999999);
            ubDict[fileInd][Ωsize[Ωl]][randNo] = [tCan,xCan,ubCan,lbCan,timedecompCan,ub5000];
            save("test_gap_ub.jld","ubDict",ubDict);
        end
    end
end


#######################################################################
# read and combine the output data
ubDict = Dict();
for fileInd in 1:4
    ubDict[fileInd] = Dict();
    for Ωl in 1:length(Ωsize)
        ubDict[fileInd][Ωsize[Ωl]] = Dict();
        for randNo in 1:20
            if randNo in keys(datagap1["ubDict"][fileInd][Ωsize[Ωl]])
                ubDict[fileInd][Ωsize[Ωl]][randNo] = datagap1["ubDict"][fileInd][Ωsize[Ωl]][randNo];
            elseif randNo in keys(datagap2["ubDict"][fileInd][Ωsize[Ωl]])
                ubDict[fileInd][Ωsize[Ωl]][randNo] = datagap2["ubDict"][fileInd][Ωsize[Ωl]][randNo];
            elseif randNo in keys(datagap3["ubDict"][fileInd][Ωsize[Ωl]])
                ubDict[fileInd][Ωsize[Ωl]][randNo] = datagap3["ubDict"][fileInd][Ωsize[Ωl]][randNo];
            elseif randNo in keys(datagap4["ubDict"][fileInd][Ωsize[Ωl]])
                ubDict[fileInd][Ωsize[Ωl]][randNo] = datagap4["ubDict"][fileInd][Ωsize[Ωl]][randNo];
            elseif randNo in keys(datagap5["ubDict"][fileInd][Ωsize[Ωl]])
                ubDict[fileInd][Ωsize[Ωl]][randNo] = datagap5["ubDict"][fileInd][Ωsize[Ωl]][randNo];
            end
        end
    end
end
# select the best solution of each test, each sample size
bestSol = Dict();
for fileInd in 1:4
    bestSol[fileInd] = Dict();
    for Ωl in 1:length(Ωsize)
        minI = indmin([ubDict[fileInd][Ωsize[Ωl]][randNo][6] for randNo in 1:20]);
        bestSol[fileInd][Ωsize[Ωl]] = [ubDict[fileInd][Ωsize[Ωl]][minI][1],ubDict[fileInd][Ωsize[Ωl]][minI][2]];
    end
end
save("test_gap_ub_sol.jld","sol",bestSol);
