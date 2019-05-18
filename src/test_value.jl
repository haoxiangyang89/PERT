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
# pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
#             "/home/haoxiang/scratch/PERT_tests/current/14/",
#             "/home/haoxiang/scratch/PERT_tests/current/19/",
#             "/home/haoxiang/scratch/PERT_tests/current/35/",
#             "/home/haoxiang/scratch/PERT_tests/current/55/",
#             "/home/haoxiang/scratch/PERT_tests/current/75/"];

# filePath = "/Users/haoxiangyang/Desktop/PERT_tests/current/14_Lognormal_Exponential/"
dDict = Dict();
ubDict = Dict();
Ωsize = 500;
global Ω = 1:Ωsize;
global ϵ = 1e-2;
global sN = 2;
global MM = 25;
global nSplit = 5;

for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize);
    global pData = pData;
    #dataRaw = load(filePath*"solData_500.jld");
    global disData = disDataSet[1];
    dDict[fileInd] = [];
    ubDict[fileInd] = [];

    global allSucc = findSuccAll(pData);
    global distanceDict = Dict();
    for i in pData.II
        for j in allSucc[i]
            distanceDict[i,j] = detCal(pData,i,j);
        end
    end

    # deterministic solution
    tic();
    tdet,xdet,fdet = detBuild(pData);
    timedet = toc();
    ubdet = ubCalP(pData,disData,Ω,xdet,tdet,999999);
    push!(dDict[fileInd],[tdet,xdet,fdet,ubdet,timedet]);
    save("test_Ext_value.jld","dDict",dDict,"ubDict",ubDict);

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
    tic();
    texp,xexp,fexp,Gexp,mexp = expModel(pData,eH,ed);
    timeexp = toc();
    ubexp = ubCalP(pData,disData,Ω,xexp,texp,999999);
    push!(dDict[fileInd],[texp,xexp,fexp,ubexp,timeexp]);
    save("test_Ext_value.jld","dDict",dDict,"ubDict",ubDict);

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
    push!(dDict[fileInd],[tdOnly,xdOnly,fdOnly,ubdOnly,timedOnly]);
    save("test_Ext_value.jld","dDict",dDict,"ubDict",ubDict);

    # HOnly solution
    for ω in Ω
        for i in pData.II
            if i != 0
                disData[ω].d[i] = mean(buildDistrn(nameD,dparams[i]));
                if disData[ω].d[i] < 1e-4
                    disData[ω].d[i] = 0;
                end
            end
        end
    end
    #tHOnly,xHOnly,fHOnly,gHOnly,mHOnly = extForm_cheat(pData,disData,Ω,1e-4,999999,noThreads);
    tic();
    #include("partSolve_Callback_tightened_sol.jl");
    tHOnly,xHOnly,ubHOnly,lbHOnly,timeIter,treeList = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,1,1e-2,6);
    timeHOnly = toc();
    gapHOnly = (ubHOnly - lbHOnly)/ubHOnly;
    fHOnly = ubHOnly;
    disData = deepcopy(disData1);
    ubHOnly = ubCalP(pData,disData,Ω,xHOnly,tHOnly,999999);
    push!(dDict[fileInd],[tHOnly,xHOnly,fHOnly,ubHOnly,gapHOnly,timeHOnly]);
    save("test_Ext_value.jld","dDict",dDict,"ubDict",ubDict);

    # full solution
    #@time tFull,xFull,fFull,gFull,mFull = extForm_cheat(pData,disData,Ω,1e-4,999999,noThreads);
    tic();
    #include("partSolve_Callback_tightened_sol.jl");
    tFull,xFull,ubFull,lbFull,timeIter,treeList = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,1,1e-2,6);
    timeFull = toc();
    gapFull = (ubCost - lbCost)/ubCost;
    ubFull = ubCost;
    fFull = lbCost;
    xFull = deepcopy(xbest);
    tFull = deepcopy(tbest);
    push!(dDict[fileInd],[tFull,xFull,fFull,ubFull,gapFull,timeFull]);
    save("test_Ext_value.jld","dDict",dDict,"ubDict",ubDict);

    dataRaw1 = load(filePath*"simuData.jld");
    Ωtest = 1:5000;
    for n in 1:20
        println("----------------Iteration $(n)----------------");
        disData1 = dataRaw1["data"][n];
        ubdet1 = ubCalP(pData,disData1,Ωtest,xdet,tdet,9999999);
        ubexp1 = ubCalP(pData,disData1,Ωtest,xexp,texp,9999999);
        ubFull1 = ubCalP(pData,disData1,Ωtest,xFull,tFull,9999999);
        ubdOnly1 = ubCalP(pData,disData1,Ωtest,xdOnly,tdOnly,9999999);
        ubHOnly1 = ubCalP(pData,disData1,Ωtest,xHOnly,tHOnly,9999999);
        push!(ubDict[fileInd],[ubdet1,ubexp1,ubFull1,ubdOnly1,ubHOnly1]);
        println(n," ",[ubdet1,ubexp1,ubFull1,ubdOnly1,ubHOnly1]);
        save("test_Ext_value.jld","dDict",dDict,"ubDict",ubDict);
    end
end
