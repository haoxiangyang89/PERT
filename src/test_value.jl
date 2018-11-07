addprocs(20);
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

pathList = ["/home/haoxiang/PERT_tests/11_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/19_Lognormal_Exponential/"];

# filePath = "/Users/haoxiangyang/Desktop/PERT_tests/14_Lognormal_Exponential/"
dDict = Dict();
for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    Ωsize = 500;
    global Ω = 1:Ωsize;
    global ϵ = 1e-2;
    pData, disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize);
    global pData = pData;
    # data = load("test_cuts.jld");
    # disData = deepcopy(data["disData"]);
    global disData = disDataSet[1];
    # dataDet = load("test_Ext_time_exponential.jld");
    allSucc = findSuccAll(pData);
    distanceDict = Dict();
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
    disData1 = deepcopy(disData);
    for ω in Ω
        disData[ω].H = mean(buildDistrn(nameH,Hparams));
    end
    tic();
    tdOnly,xdOnly,fdOnly,gdOnly,mdOnly = extForm_cheat(pData,disData,Ω,1e-4,999999);
    timedOnly = toc();
    disData = deepcopy(disData1);
    ubdOnly = ubCal(pData,disData,Ω,xdOnly,tdOnly,999999);

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
    tic();
    include("partSolve_Callback_tightened.jl");
    timeHOnly = toc();
    gapHOnly = (ubCost - lbCost)/ubCost;
    xHOnly = deepcopy(xbest);
    tHOnly = deepcopy(tbest);
    disData = deepcopy(disData1);
    ubHOnly = ubCal(pData,disData,Ω,xHOnly,tHOnly,999999);

    # full solution
    tic();
    include("partSolve_Callback_tightened.jl");
    timeFull = toc();
    gapFull = (ubCost - lbCost)/ubCost;
    ubFull = ubCost;
    lbFull = lbCost;
    xFull = deepcopy(xbest);
    tFull = deepcopy(tbest);

    dDict[fileInd] = [tdet,xdet,fdet,ubdet,
                texp,xexp,fexp,Gexp,ubexp,
                tFull,xFull,ubFull,lbFull,
                tdOnly,xdOnly,gdOnly,ubdOnly,
                tHOnly,xHOnly,ubHOnly,gapHOnly];

    ubList = [];
    for n in 1:30
        println("----------------Iteration $(n)----------------");
        disData1,Ω = autoUGen("LogNormal",Hparams,"Exponential",dparams,2000,1 - pData.p0);
        disData1 = orderdisData(disData1,Ω);
        ubdet1 = ubCal(pData,disData1,Ω,xdet,tdet,999999);
        ubexp1 = ubCal(pData,disData1,Ω,xexp,texp,999999);
        ubFull1 = ubCal(pData,disData1,Ω,xFull,tFull,999999);
        ubdOnly1 = ubCal(pData,disData1,Ω,xdOnly,tdOnly,999999);
        ubHOnly1 = ubCal(pData,disData1,Ω,xHOnly,tHOnly,999999);
        push!(ubList,[ubdet1,ubexp1,ubFull1,ubdOnly1,ubHOnly1]);
        println(n," ",[ubdet1,ubexp1,ubFull1,ubdOnly1,ubHOnly1]);
    end
    save("test_Ext_value.jld","dDict",dDict,"ubList",ubList);
end
