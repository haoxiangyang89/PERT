addprocs(20);
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

filePath = "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/";
Ωsize = [10,50,100,200,500,1000];
dDict = Dict();
for Ωl in 1:length(Ωsize)
    Ω = 1:Ωsize[Ωl];
    ϵ = 1e-2;
    pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize[Ωl]);
    disData = disDataSet[1];

    allSucc = findSuccAll(pData);
    distanceDict = Dict();
    for i in pData.II
        for j in allSucc[i]
            distanceDict[i,j] = detCal(pData,i,j);
        end
    end
    # our decomposition method
    tic();
    include("partSolve_Callback_tightened.jl");
    timedecomp = toc();
    gapdecomp = (ubCost - lbCost)/ubCost;

    dDict[Ωsize[Ωl]] = [tbest,xbest,lbCost,ubCost,gapdecomp,timedecomp];
end
save("test_Ext_budget.jld","dDict",dDict);

################################################################
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

filePath = "/home/haoxiang/PERT_tests/14_Lognormal_Exponential";
Ωsize = 500;
Ω = 1:Ωsize;
ϵ = 1e-2;
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize);
data_value = load("test_Ext_value.jld");
tdet = data_value["dDict"][1];
xdet = data_value["dDict"][2];

texp = data_value["dDict"][5];
xexp = data_value["dDict"][6];

tFull = data_value["dDict"][10];
xFull = data_value["dDict"][11];

tdOnly = data_value["dDict"][14];
xdOnly = data_value["dDict"][15];

tHOnly = data_value["dDict"][18];
xHOnly = data_value["dDict"][19];

ubList = [];
for n in 1:30
    println("----------------Iteration $(n)----------------");
    disData1,Ω = autoUGen("LogNormal",Hparams,"Exponential",dparams,500,1 - pData.p0);
    disData1 = orderdisData(disData1,Ω);
    ubdet1 = ubCal(pData,disData1,Ω,xdet,tdet,999999);
    ubexp1 = ubCal(pData,disData1,Ω,xexp,texp,999999);
    ubFull1 = ubCal(pData,disData1,Ω,xFull,tFull,999999);
    ubdOnly1 = ubCal(pData,disData1,Ω,xdOnly,tdOnly,999999);
    ubHOnly1 = ubCal(pData,disData1,Ω,xHOnly,tHOnly,999999);
    push!(ubList,[ubdet1,ubexp1,ubFull1,ubdOnly1,ubHOnly1]);
    println(n," ",[ubdet1,ubexp1,ubFull1,ubdOnly1,ubHOnly1]);
end
save("test_Ext_value_Out.jld","ubList",ubList);
