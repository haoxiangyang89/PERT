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
Ωsize = [10,50,100,200,500,1000];
for fileInd in 1:length(pathList)
    dDict[fileInd] = Dict();
    filePath = pathList[fileInd];
    pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,5000);
    global Ω1 = 1:5000;
    global disData1 = disDataSet[1];
    for Ωl in 1:length(Ωsize)
        global Ω = 1:Ωsize[Ωl];
        global ϵ = 1e-2;
        dDict[fileInd][Ωsize[Ωl]] = [];
        ubbudget = [];
        n = 1;
        while n <= 20
            try
            pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize[Ωl]);
            global disData = disDataSet[1];
            global pData = pData;

            global allSucc = findSuccAll(pData);
            global distanceDict = Dict();
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
            xFull = deepcopy(xbest);
            tFull = deepcopy(tbest);

            ubTemp = ubCalP(pData,disData1,Ω1,xFull,tFull,999999);
            push!(dDict[fileInd][Ωsize[Ωl]],[tFull,xFull,lbCost,ubCost,gapdecomp,timedecomp,ubTemp]);
            save("test_Ext_budget.jld","dDict",dDict);
            n += 1;
            catch
            println("Error in Data!");
            end
        end
    end
end
