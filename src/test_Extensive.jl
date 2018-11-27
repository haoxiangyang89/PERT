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
        # extensive formulation
        tic();
        text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,ϵ,999999);
        timeext = toc();
        ubmp = mext.objVal;
        lbmp = mext.objBound;
        gapext = (mext.objVal - mext.objBound)/mext.objVal;
        ubext = ubCal(pData,disData,Ω,xext,text,999999);
        dDict[fileInd][Ωsize[Ωl]] = [text,xext,lbmp,ubmp,ubext,gapext,timeext];
        save("test_Ext_extensive.jld","dDict",dDict);
    end
end
