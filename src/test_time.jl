addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

#filePath = "/Users/haoxiangyang/Desktop/PERT_tests/75_Lognormal_Exponential/"
pathList = ["/home/haoxiang/PERT_tests/11_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/19_Lognormal_Exponential/"];
Ωsize = [100,200,300,400,500,750,1000];
sNList = [20,20,20,20,20,30,40];
MMList = [5,10,15,20,25,25,25];
dDict = Dict();
for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    dDict[fileInd] = Dict();
    for Ωl in 1:length(Ωsize)
        global Ω = 1:Ωsize[Ωl];
        global ϵ = 1e-2;
        global pData;
        pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize[Ωl]);
        global disData = disDataSet[1];

        global allSucc = findSuccAll(pData);
        global distanceDict = Dict();
        for i in pData.II
            for j in allSucc[i]
                distanceDict[i,j] = detCal(pData,i,j);
            end
        end
        # our decomposition method
        global sN = sNList[Ωl];
        global MM = MMList[Ωl];
        tic();
        include("partSolve_Callback_tightened_sol.jl");
        timedecomp = toc();
        gapdecomp = (ubCost - lbCost)/ubCost;
        ubFull = ubCost;
        lbFull = lbCost;
        xFull = deepcopy(xbest);
        tFull = deepcopy(tbest);

        tic();
        include("partSolve_Callback_tightened.jl");
        timedecomp1 = toc();
        gapdecomp1 = (ubCost - lbCost)/ubCost;
        ubFull1 = ubCost;
        lbFull1 = lbCost;
        xFull1 = deepcopy(xbest);
        tFull1 = deepcopy(tbest);

        # extensive formulation
        tic();
        text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,1e-2,999999);
        timeext = toc();
        ubmp = mext.objVal;
        lbmp = mext.objBound;
        gapext = (mext.objVal - mext.objBound)/mext.objVal;
        ubext = ubCal(pData,disData,Ω,xext,text,999999);
        dDict[fileInd][Ωsize[Ωl]] = [tFull,xFull,lbFull,ubFull,gapdecomp,timedecomp,
                            tFull1,xFull1,lbFull1,ubFull1,gapdecomp1,timedecomp1,
                            text,xext,lbmp,ubmp,ubext,gapext,timeext];
        save("test_Ext_time.jld","dDict",dDict);
    end
end
