# output the result
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

datalb1 = load("test_gap_lb_1.jld");
datalb2 = load("test_gap_lb_2.jld");
datalb3 = load("test_gap_lb_3.jld");
datalb4 = load("test_gap_lb_4.jld");
datalb5 = load("test_gap_lb_5.jld");

lbDict = Dict();
for fileInd in 1:4
    lbDict[fileInd] = Dict();
    for Ωl in 1:length(Ωsize)
        lbDict[fileInd][Ωsize[Ωl]] = [];
        append!(lbDict[fileInd][Ωsize[Ωl]],datalb1["lbDict"][fileInd][Ωsize[Ωl]]);
        append!(lbDict[fileInd][Ωsize[Ωl]],datalb2["lbDict"][fileInd][Ωsize[Ωl]]);
        append!(lbDict[fileInd][Ωsize[Ωl]],datalb3["lbDict"][fileInd][Ωsize[Ωl]]);
        append!(lbDict[fileInd][Ωsize[Ωl]],datalb4["lbDict"][fileInd][Ωsize[Ωl]]);
        append!(lbDict[fileInd][Ωsize[Ωl]],datalb5["lbDict"][fileInd][Ωsize[Ωl]]);
    end
end

outputData = Dict();
for fileInd in 1:4
    filePath = pathList[fileInd]
    simuPath = filePath*"simuData.jld";
    data5000Raw = load(simuPath);
    outputData[fileInd] = Dict();

    for Ωl in 1:length(Ωsize)
        tCanBest = bestSol[fileInd][Ωsize[Ωl]][1];
        xCanBest = bestSol[fileInd][Ωsize[Ωl]][2];
        # obtain U5000
        ubList = [];
        for randNo in 1:20
            disData1 = data5000Raw["data"][randNo];
            Ω1 = 1:length(disData1);
            ub5000 = ubCalP(pData,disData1,Ω1,xCanBest,tCanBest,999999);
            push!(ubList,ub5000);
        end

        gapList = [lbDict[fileInd][Ωsize[Ωl]][randNo][6] - lbDict[fileInd][Ωsize[Ωl]][randNo][3] for randNo in 1:20]
        gapmean = mean(gapList);
        gapsigma = std(gapList);
        gapub = gapmean + 1.729*gapsigma/sqrt(20);
        ubmean = mean(ubList);
        gapratio = gapub/ubmean*100;
        outputData[fileInd][Ωsize[Ωl]] = [gapmean,gapub,gapratio];
    end
end
save("test_gap_out.jld","data",outputData);
