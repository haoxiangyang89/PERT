addprocs(20);
# test sbb
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("master.jl");
@everywhere include("sub.jl");
@everywhere include("cuts.jl");
@everywhere include("iSolve.jl");
@everywhere include("tighten.jl");
@everywhere include("branchFunc.jl");
@everywhere include("detForm.jl");
@everywhere include("extForm.jl");
@everywhere include("ubCalFunc.jl");
@everywhere include("tighten.jl");
@everywhere include("partition_LP.jl");
@everywhere include("partition_LR.jl");
@everywhere include("part_tight.jl");
@everywhere include("partitionSolve.jl");

pInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_P.csv";
kInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_K.csv";
ϕInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_Phi_full.csv";

pData = readInP(pInputAdd,kInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
disData,Ω = autoUGen("LogNormal",[log(35),0.5],nameD,dparams,500,1 - pData.p0);
disData = orderdisData(disData,Ω);

tic();
tbest,xbest,lbCost,ubCost = partitionSolve(pData,disData,0.01);
timedecomp = toc();
gapdecomp = (ubCost - lbCost)/ubCost;

tic();
text,xext,fext,gext,mp = extForm_cheat(pData,disData,Ω,1e-2,999999);
timeext = toc();
ubmp = mp.objVal;
lbmp = mp.objBound;
gapext = (mp.objVal - mp.objBound)/mp.objVal;
dDict = [tbest,xbest,lbCost,ubCost,gapdecomp,timedecomp,text,xext,fext,gext,ubmp,lbmp,gapext,timeext];
save("test_Ext_time.jld","dDict",dDict);
