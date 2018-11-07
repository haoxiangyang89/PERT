addprocs(20);
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

#filePath = "/Users/haoxiangyang/Desktop/PERT_tests/75_Lognormal_Exponential/"
pathList = ["/home/haoxiang/PERT_tests/11_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/19_Lognormal_Exponential/"];

for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    Ωsize = [100,200,500,1000];
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
        global sN = 20;
        global MM = 25;
        tic();
        include("partSolve_Callback_tightened_sol.jl");
        timedecomp = toc();
        gapdecomp = (ubCost - lbCost)/ubCost;

        # extensive formulation
        tic();
        text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,1e-2,999999);
        timeext = toc();
        ubmp = mext.objVal;
        lbmp = mext.objBound;
        gapext = (mext.objVal - mext.objBound)/mext.objVal;
        ubext = ubCal(pData,disData,Ω,xext,text,999999);
        dDict[Ωsize[Ωl]] = [tbest,xbest,lbCost,ubCost,gapdecomp,timedecomp,
                            text,xext,lbmp,ubmp,gapext,timeext];
    end
    save("test_Ext_time.jld","dDict",dDict);
end
