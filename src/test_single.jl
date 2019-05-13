# test the running time given a number of scenarios and obtain the result

addprocs(3);
global noThreads = 4;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

# pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
#             "/home/haoxiang/scratch/PERT_tests/current/14/",
#             "/home/haoxiang/scratch/PERT_tests/current/19/",
#             "/home/haoxiang/scratch/PERT_tests/current/35/",
#             "/home/haoxiang/scratch/PERT_tests/current/55/",
#             "/home/haoxiang/scratch/PERT_tests/current/75/"];

# pathList = ["/home/haoxiang/PERT_tests/current/11/",
#             "/home/haoxiang/PERT_tests/current/14/",
#             "/home/haoxiang/PERT_tests/current/19/",
#             "/home/haoxiang/PERT_tests/current/35/",
#             "/home/haoxiang/PERT_tests/current/55/",
#             "/home/haoxiang/PERT_tests/current/75/"];

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
ΩsizeTest = 500;
global Ωt = 1:ΩsizeTest;
global ϵ = 1e-4;
global pData;
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize,20);
pData,disDataSetT,nameD,nameH,dparams,Hparams = genData(filePath,ΩsizeTest,20);

dataList = [];
for itemI in 1:20
    global pData = pData;
    global disData = disDataSet[itemI];
    global disDataTest = disDataSetT[itemI];

    global allSucc = findSuccAll(pData);
    global distanceDict = Dict();
    for i in pData.II
        for j in allSucc[i]
            distanceDict[i,j] = detCal(pData,i,j);
        end
    end

    tdet,xdet,fdet = detBuild(pData);
    ubdet = ubCalP(pData,disDataTest,Ωt,xdet,tdet,999999);

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
    ubexp,ubexpList = ubCalP(pData,disDataTest,Ωt,xexp,texp,999999,1);

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
    ubdOnly,ubdOList = ubCalP(pData,disDataTest,Ωt,xdOnly,tdOnly,999999,1);

    tic();
    text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,1e-4,999999,noThreads);
    timeext = toc();
    ubmp,ubmpList = ubCalP(pData,disDataTest,Ωt,xext,text,999999,1);

    push!(dataList,[ubdet,ubexp,ubdOnly,ubmp]);
end

save("test_single.jld","data",[text,xext,fext,gext]);
