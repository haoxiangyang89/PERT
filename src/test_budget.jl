addprocs(20);
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

filePath = "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/";
Ωsize = [10,50,100,200,500,1000];
dDict = Dict();
for Ωl in 5:length(Ωsize)
    Ω = 1:Ωsize[Ωl];
    ϵ = 1e-2;
    dDict[Ωsize[Ωl]] = [];
    n = 1;
    while n <= 20
        try
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

            push!(dDict[Ωsize[Ωl]],[tbest,xbest,lbCost,ubCost,gapdecomp,timedecomp]);
            save("test_Ext_budget.jld","dDict",dDict);
            n += 1;
        catch
            println("Error in Solving Process!");
        end
    end
end

################################################################
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

filePath = "/home/haoxiang/PERT_tests/14_Lognormal_Exponential";
ϵ = 1e-2;
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,500);
data_budget = load("test_Ext_budget.jld");
xFull = Dict();
tFull = Dict();
ubbudget = Dict();
for n in 1:20
    xFull[n] = Dict();
    tFull[n] = Dict();
    ubbudget[n] = Dict();
    for Ωl in 1:length(Ωsize)
        if n in keys(data_budget["dDict"][Ωsize[Ωl]])
            xFull[n][Ωsize[Ωl]] = data_budget["dDict"][Ωsize[Ωl]][n][2];
            tFull[n][Ωsize[Ωl]] = data_budget["dDict"][Ωsize[Ωl]][n][1];
        end
    end
end

disData1,Ω = autoUGen("LogNormal",Hparams,"Exponential",dparams,5000,1 - pData.p0);
disData1 = orderdisData(disData1,Ω);
ubTempList = [];
for Ωl in 1:length(Ωsize)
    for n in 1:20
        if Ωsize[Ωl] in keys(xFull[n])
            ubTemp = ubCal(pData,disData1,Ω,xFull[n][Ωsize[Ωl]],tFull[n][Ωsize[Ωl]],999999);
            ubbudget[n][Ωsize[Ωl]] = ubTemp;
        end
    end
end
save("test_Ext_budget_Out.jld","ubbudget",ubbudget);
