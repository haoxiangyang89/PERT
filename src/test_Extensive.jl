using Distributed;
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,Ipopt,MathProgBase;
@everywhere using Distributions,HDF5,JLD,DelimitedFiles,Statistics,SharedArrays;
# test sbb
@everywhere include("header.jl");

#filePath = "/Users/haoxiangyang/Desktop/PERT_tests/75_Lognormal_Exponential/"
pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
            "/home/haoxiang/scratch/PERT_tests/current/14/",
            "/home/haoxiang/scratch/PERT_tests/current/19/",
            "/home/haoxiang/scratch/PERT_tests/current/35/",
            "/home/haoxiang/scratch/PERT_tests/current/55/",
            "/home/haoxiang/scratch/PERT_tests/current/75/"];
# pathList = ["/scratch/haoxiang/current/11/",
#                "/scratch/haoxiang/current/14/",
#                "/scratch/haoxiang/current/19/",
#                "/scratch/haoxiang/current/35/",
#                "/scratch/haoxiang/current/55/",
#                "/scratch/haoxiang/current/75/"];

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
        tempTimer = time();
        text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,ϵ,999999);
        timeext = time() - tempTimer;
        ubmp = mext.objVal;
        lbmp = mext.objBound;
        gapext = (mext.objVal - mext.objBound)/mext.objVal;
        ubext = ubCal(pData,disData,Ω,xext,text,999999);
        dDict[fileInd][Ωsize[Ωl]] = [text,xext,lbmp,ubmp,ubext,gapext,timeext];
        save("test_Ext_extensive.jld","dDict",dDict);
    end
end
