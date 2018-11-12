addprocs(30);
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

#filePath = "/Users/haoxiangyang/Desktop/PERT_tests/75_Lognormal_Exponential/"
pathList = ["/home/haoxiang/PERT_tests/11_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/19_Lognormal_Exponential/"];

ErrorData = [];
dDict = Dict();
for fileInd in 1:length(pathList)
    dDict[fileInd] = Dict();
    filePath = pathList[fileInd];
    Ωsize = [100,200,300,400,500,750,1000];
    for Ωl in 1:length(Ωsize)
        global Ω = 1:Ωsize[Ωl];
        global ϵ = 1e-2;
        global pData;
        pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize[Ωl]);
        global pData = pData;
        global disData = disDataSet[1];

        global allSucc = findSuccAll(pData);
        global distanceDict = Dict();
        for i in pData.II
            for j in allSucc[i]
                distanceDict[i,j] = detCal(pData,i,j);
            end
        end
        # our decomposition method
        tic();
        try
            include("partSolve_Callback_tightened.jl");
        catch
            push!(ErrorData,disData);
            save("ErrorData.jld","data",ErrorData);
        end
        timedecomp = toc();
        gapdecomp = (ubCost - lbCost)/ubCost;
        ubFull = ubCost;
        lbFull = lbCost;
        xFull = deepcopy(xbest);
        tFull = deepcopy(tbest);

        # extensive formulation
        tic();
        text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,1e-2,999999);
        timeext = toc();
        ubmp = mext.objVal;
        lbmp = mext.objBound;
        gapext = (mext.objVal - mext.objBound)/mext.objVal;
        ubext = ubCal(pData,disData,Ω,xext,text,999999);
        dDict[fileInd][Ωsize[Ωl]] = [tFull,xFull,lbFull,ubFull,gapdecomp,timedecomp,
                            text,xext,lbmp,ubmp,ubext,gapext,timeext];
        save("test_Ext_time.jld","dDict",dDict);
    end
end
