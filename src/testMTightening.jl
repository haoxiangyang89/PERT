using Distributed;
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,Ipopt,MathProgBase;
@everywhere using Distributions,HDF5,JLD,DelimitedFiles,Statistics,SharedArrays;
@everywhere const GUROBI_ENV = Gurobi.Env();

# test sbb
@everywhere include("header.jl");

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

Ωsize = [100,200,500,1000];
dDict = Dict();
ϵ = 1e-2;
sNList = [10,20,25,25];
MMList = [10,10,20,40];

fileInd = 1;
filePath = pathList[fileInd];
Ωl = 1;
global Ω = 1:Ωsize[Ωl];
randNo = 1;
extBool = true;

pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1);
global pData = pData;
disDataRaw = load(pathList[fileInd]*"solData_100.jld");
disData = disDataRaw["data"][randNo];

global allSucc = findSuccAll(pData);
global distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

global sN = 20;
global MM = 5;
text,xext,fext,gext,mext = extForm_cheat_new(pData,disData,Ω,sN,MM,1e-2,99999,noThreads);

#%%
# Case 14 sample size 200
fileInd = 2;
filePath = pathList[fileInd];
Ωl = 2;
global Ω = 1:Ωsize[Ωl];
randNo = 1;
extBool = true;

pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1);
global pData = pData;
disDataRaw = load(pathList[fileInd]*"solData_$(Ωsize[Ωl]).jld");
disData = disDataRaw["data"][randNo];

global allSucc = findSuccAll(pData);
global distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end

global sN = sNList[Ωl];
global MM = MMList[Ωl];
text,xext,fext,gext,mext,timeext = extForm_cheat_new(pData,disData,Ω,sN,MM,1e-2,99999,noThreads);

Mω = getMomega(pData,disData);
text_u,xext_u,fext_u,gext_u,mext_u,timeext_u = extForm_cheat_new(pData,disData,Ω,sN,MM,1e-2,99999,noThreads,"",Mω);
