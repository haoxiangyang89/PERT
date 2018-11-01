# create data set for each of the problem and store them as .jld files
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

function createData(filePath,Ωsize,jldName = "testData_$(Ωsize).jld")
    jldPath = joinpath(filePath,jldName);
    Ω = 1:Ωsize;
    pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize);
    save(jldPath,"disDataSet",disDataSet);
end

pathList = ["/home/haoxiang/PERT_tests/11_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/19_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/55_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/75_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/110_Lognormal_Exponential/"];

for path in pathList
    createData(path,500);
end

# test data quality
sN = 20;
MM = 25;
ubList,tHList = iniPart(pData,disData,Ω,sN,MM);
