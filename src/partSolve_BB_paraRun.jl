addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");
@everywhere include("partSolve_BB_shared.jl");

pathList = ["/home/haoxiang/PERT_tests/11_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/19_Lognormal_Exponential/"];
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
# data = load("test_cuts.jld");
# disData = deepcopy(data["disData"]);
global disData = disDataSet[1];
# dataDet = load("test_Ext_time_exponential.jld");
global allSucc = findSuccAll(pData);
global distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

# deterministic solution
tdet,xdet,fdet = detBuild(pData);
ubdet = ubCalP(pData,disData,Ω,xdet,tdet,999999);

# expected solution
eH = mean(buildDistrn(nameH,Hparams));
ed = Dict();
for i in pData.II
    if i == 0
        ed[i] = 0;
    else
        ed[i] = mean(buildDistrn(nameD,dparams[i]));
    end
end
texp,xexp,fexp,Gexp,mexp = expModel(pData,eH,ed);
ubexp = ubCalP(pData,disData,Ω,xexp,texp,999999);

# dOnly solution
global disData1;
disData1 = deepcopy(disData);
for ω in Ω
    disData[ω].H = mean(buildDistrn(nameH,Hparams));
end
tic();
tdOnly,xdOnly,fdOnly,gdOnly,mdOnly = extForm_cheat(pData,disData,Ω,1e-4,999999,noThreads);
timedOnly = toc();
disData = deepcopy(disData1);
ubdOnly = ubCalP(pData,disData,Ω,xdOnly,tdOnly,999999);

global sN = 25;
global MM = 20;
global r = 1e-6;
tic();
data = partSolve_BB_para(pData,disData,Ω,sN,MM,noThreads,1e-4);
timeLast = toc();
save("testBB_data.jld","data",data,"timeLast",timeLast);
