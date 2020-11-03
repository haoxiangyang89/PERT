using Distributed;
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,Ipopt,MathProgBase;
@everywhere using Distributions,HDF5,JLD,DelimitedFiles,Statistics,SharedArrays;
@everywhere const GUROBI_ENV = Gurobi.Env();

# test sbb
@everywhere include("header.jl");
@everywhere include("pullDecomp.jl");

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
Ωsize = [200];
Ωl = 5;
global Ω = 1:Ωsize[Ωl];

dDict = Dict();
for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    bestObj,that,xhat = pullDecomp(pData,disData,Ω,1e-2);
    dDict[fileInd] = [bestObj,that,xhat];
    save("test_pullDecomp.jld","data",dDict);
end
