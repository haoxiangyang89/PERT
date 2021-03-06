# test the extensive formulation with a preset upper bound
using Distributed;
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,Ipopt,MathProgBase;
@everywhere using Distributions,HDF5,JLD,DelimitedFiles,Statistics,SharedArrays;
@everywhere const GUROBI_ENV = Gurobi.Env()

# test sbb
@everywhere include("header.jl");

pathList = ["/scratch/haoxiang/current/11/",
            "/scratch/haoxiang/current/14/",
            "/scratch/haoxiang/current/19/",
            "/scratch/haoxiang/current/35/",
            "/scratch/haoxiang/current/55/",
            "/scratch/haoxiang/current/75/"];
# pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
#             "/home/haoxiang/scratch/PERT_tests/current/14/",
#             "/home/haoxiang/scratch/PERT_tests/current/19/",
#             "/home/haoxiang/scratch/PERT_tests/current/35/",
#             "/home/haoxiang/scratch/PERT_tests/current/55/",
#             "/home/haoxiang/scratch/PERT_tests/current/75/"];

# pathList = ["/Users/haoxiangyang/Desktop/PERT_tests/current/11_Lognormal_Exponential",
#             "/Users/haoxiangyang/Desktop/PERT_tests/current/14_Lognormal_Exponential",
#             "/Users/haoxiangyang/Desktop/PERT_tests/current/19_Lognormal_Exponential"]
Ωsize = [100,200,500,1000];
sNList = [20,20,25,25];
MMList = [5,10,20,40];
data = Dict();

for fileInd in 1:4
    filePath = pathList[fileInd];
    data[fileInd] = Dict();

    pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1);
    global pData = pData;
    global allSucc = findSuccAll(pData);
    global distanceDict = Dict();
    for i in pData.II
        for j in allSucc[i]
            distanceDict[i,j] = detCal(pData,i,j);
        end
    end
    for Ωl in 1:4
        global Ω = 1:Ωsize[Ωl];
        disDataRaw = load(pathList[fileInd]*"solData_$(Ωsize[Ωl]).jld");
        global sN = sNList[Ωl];
        global MM = MMList[Ωl];
        dataList = [];
        for randNo in 1:20
            disData = disDataRaw["data"][randNo];
            tempTimer = time();
            ubextList,tHList,ubInc,tbest,xbest,θbest,textList,xextList = iniPart(pData,disData,Ω,sN,MM,1,noThreads);
            tub = time() - tempTimer;
            push!(dataList,(ubInc,tub));
        end
        data[fileInd][Ωsize[Ωl]] = dataList;
        save("test_UB.jld","data",data);
    end
end

################################################################################
# print the output
dataB = load("test_Ext_budget.jld");
gapDict = Dict();
timeDict = Dict();
for i in 1:4
    for j in [100,200,500,1000]
        gapDict[i,j] = [(data["data"][i][j][k][1] - dataB["dDict"][i][j][k][4])/dataB["dDict"][i][j][k][4]*100 for k in 1:20];
        timeDict[i,j] = [data["data"][i][j][k][2] for k in 1:20];
    end
end
