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
# pathList = ["/Users/haoxiangyang/Desktop/PERT_tests/current/11_Lognormal_Exponential",
#             "/Users/haoxiangyang/Desktop/PERT_tests/current/14_Lognormal_Exponential",
#             "/Users/haoxiangyang/Desktop/PERT_tests/current/19_Lognormal_Exponential"]
fileInd = 2;
filePath = pathList[fileInd];
Ωsize = 500;
global Ω = 1:Ωsize;
global ϵ = 1e-4;
global pData;
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize);
global pData = pData;
dataRaw = load("14_test1.jld");
# disData = deepcopy(data["disData"]);
global disData = dataRaw["disData"];
# dataDet = load("test_Ext_time_exponential.jld");
global allSucc = findSuccAll(pData);
global distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

global sN = 25;
global MM = 20;
global r = 1e-6;
tic();
data = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,1e-4);
timeLast = toc();
save("testBB_data.jld","data",data,"timeLast",timeLast);
