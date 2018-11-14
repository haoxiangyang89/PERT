addprocs(30);
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

pathList = ["/home/haoxiang/PERT_tests/11_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/19_Lognormal_Exponential/"];

# filePath = "/Users/haoxiangyang/Desktop/PERT_tests/14_Lognormal_Exponential/"
dDict = Dict();
Ωsize = [10,20,50,75,100,200,300,400,500];
sNList = [0,0,0,0,20,20,20,20,20];
MMList = [0,0,0,0,5,10,15,20,25];
data5000Raw = load("data5000.jld");
data1 = data5000Raw["dataUB"];
for fileInd in 1:length(pathList)
    dDict[fileInd] = Dict();
    filePath = pathList[fileInd];
    disData1 = data1[fileInd];
    for Ωl in 1:length(Ωsize)
        global Ω = 1:Ωsize1[Ωl];
        global ϵ = 1e-2;
        dDict[fileInd][Ωsize1[Ωl]] = [];
        ubbudget = [];
        n = 1;
        while n <= 20
            # try
            global pData;
            global disDataSet;
            pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize1[Ωl]);
            global disData = disDataSet[1];

            global allSucc = findSuccAll(pData);
            global distanceDict = Dict();
            for i in pData.II
                for j in allSucc[i]
                    distanceDict[i,j] = detCal(pData,i,j);
                end
            end
            # our decomposition method
            if Ωl < 5
                tic();
                include("partSolve_Callback_tightened.jl");
                timedecomp = toc();
                gapdecomp = (ubCost - lbCost)/ubCost;
                xFull = deepcopy(xbest);
                tFull = deepcopy(tbest);
            else
                global sN = sNList[Ωl];
                global MM = MMList[Ωl];
                tic();
                include("partSolve_Callback_tightened_sol.jl");
                timedecomp = toc();
                gapdecomp = (ubCost - lbCost)/ubCost;
                xFull = deepcopy(xbest);
                tFull = deepcopy(tbest);
            end

            ubTemp = ubCalP(pData1,disData1,Ω1,xFull,tFull,999999);
            push!(dDict[fileInd][Ωsize1[Ωl]],[tFull,xFull,lbCost,ubCost,gapdecomp,timedecomp,ubTemp]);
            save("test_Ext_budget.jld","dDict",dDict);
            n += 1;
            # catch
            # println("Error in Data!");
            # push!(ErrorData,(fileInd,disData));
            # end
        end
    end
end
