# script to make tests
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
@everywhere include("header.jl");

pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
            "/home/haoxiang/scratch/PERT_tests/current/14/",
            "/home/haoxiang/scratch/PERT_tests/current/19/",
            "/home/haoxiang/scratch/PERT_tests/current/35/",
            "/home/haoxiang/scratch/PERT_tests/current/55/",
            "/home/haoxiang/scratch/PERT_tests/current/75/"];

ΩsizeSet = [10,20,50,100,200,300,400,500,750,1000,1500,2000];
for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    for Ωsize in ΩsizeSet
        Ω = 1:Ωsize;
        ϵ = 1e-4;
        pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize,20);
        save(filePath*"solData_$Ωsize.jld","data",disDataSet);
    end
    ΩsizeTest = 5000;
    Ωt = 1:ΩsizeTest;
    pData,disDataSetT,nameD,nameH,dparams,Hparams = genData(filePath,ΩsizeTest,20);
    save(filePath*"simuData.jld","data",disDataSetT);
end
