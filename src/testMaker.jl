# script to make tests
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
@everywhere include("header.jl");

pathList = ["/home/haoxiang/PERT_tests/current/11/",
            "/home/haoxiang/PERT_tests/current/14/",
            "/home/haoxiang/PERT_tests/current/19/",
            "/home/haoxiang/PERT_tests/current/35/",
            "/home/haoxiang/PERT_tests/current/55/",
            "/home/haoxiang/PERT_tests/current/75/"];

# make the skewness test
Dmag = 1000;
μp = 1/2;
HsList = 0.1:0.1:0.9;
dataSize = 50;
Ωsize = 100;
for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    for Hskewness in HsList
        pData,disDataSet = experimentgen(filePath,dataSize,Ωsize,Hskewness,μp,Dmag);
        save(filePath*"test1_$(fileInd)_$(Hskewness).jld","pData",pData,"disDataSet",disDataSet);
    end
end

# make the d test
Dmag = 1000;
μp = 1/2;
Hskewness = 1/2;
dataSize = 50;
Ωsize = 100;
for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    for dOpt in ["a","d","r"]
        pData,disDataSet = experimentgen(filePath,dataSize,Ωsize,Hskewness,μp,Dmag,dOpt);
        save(filePath*"test2_$(fileInd)_$(dOpt).jld","pData",pData,"disDataSet",disDataSet);
    end
end
