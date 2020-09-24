# test whether the box formed by solve small sample problem will contain the optimal solution
using Distributed;
addprocs(3);
global noThreads = 4;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

pathList = ["/Users/haoxiangyang/Desktop/PERT_tests/current/11/",
            "/Users/haoxiangyang/Desktop/PERT_tests/current/14/",
            "/Users/haoxiangyang/Desktop/PERT_tests/current/19/",
            "/Users/haoxiangyang/Desktop/PERT_tests/current/35/",
            "/Users/haoxiangyang/Desktop/PERT_tests/current/55/",
            "/Users/haoxiangyang/Desktop/PERT_tests/current/75/"];

fileInd = 2;
filePath = pathList[fileInd];
Ωsize = 50;
global Ω = 1:Ωsize;
global ϵ = 1e-4;
global pData;
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize,20,"test_P.csv","test_K.csv","test_Phi_unif.csv");

# start with an upper bound based on the smaller stochastic solution
sN = 2;
MM = 25;
dataList = [];

for iterN in 1:20
    disData = disDataSet[iterN];

    H = Dict();
    H[0] = 0;
    # remove the duplicate ones
    counter = 1;
    for ω in Ω
        if !(disData[ω].H in values(H))
            H[counter] = disData[ω].H;
            counter += 1;
        end
    end
    Tmax = disData[length(Ω)].H + longestPath(pData)[0];
    H[counter] = Tmax;

    ubextList,tHList,ubInc,tbest,xbest,θbest,textList,xextList = iniPart(pData,disData,Ω,sN,MM,1,noThreads);
    text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,1e-4,999999,noThreads);
    for item in tHList
        if (text[item[1]] > H[item[3]])|(text[item[1]] < H[item[2]])
            push!(dataList,[iterN,tHList,text,xext,fext,gext]);
            println(item[1]," ",text[item[1]]," ",H[item[2]]," ",H[item[3]]);
        end
    end
end
