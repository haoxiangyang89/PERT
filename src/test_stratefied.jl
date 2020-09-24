using Distributed;
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

pathList = ["/home/haoxiang/PERT_tests/11_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/14_Lognormal_Exponential/",
            "/home/haoxiang/PERT_tests/19_Lognormal_Exponential/"];

strataList = [250,100,50,25,10,5,2];
sampleList = [2,5,10,20,50,100,250];
data5000Raw = load("data5000.jld");
data1 = data5000Raw["dataUB"];
dDict = Dict();
for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    global ϵ = 1e-2;
    global pData;
    dDict[fileInd] = Dict();
    disData1 = data1[fileInd];
    Ω1 = 1:length(disData1);
    for j in 1:length(strataList)
        dDict[fileInd][strataList[j],sampleList[j]] = [];
        for n in 1:20
            pData, disDataSet,nameD,nameH,dparams,Hparams = genDataStrata(filePath,strataList[j],sampleList[j]);
            global pData = pData;
            global disData = disDataSet[1];
            Ω = 1:(strataList[j]*sampleList[j]);
            tic();
            text,xext,fext,gext,mext = extForm_cheat_new(pData,disData,Ω,1e-2,999999);
            timeext = toc();
            ubmp = mext.objVal;
            lbmp = mext.objBound;
            gapext = (mext.objVal - mext.objBound)/mext.objVal;
            ubTemp = ubCalP(pData,disData1,Ω1,xext,text,999999);
            push!(dDict[fileInd][strataList[j],sampleList[j]],[text,xext,fext,gext,ubmp,lbmp,gapext,ubTemp]);
            save("test_Ext_stratefied.jld","dDict",dDict);
        end
    end
end
