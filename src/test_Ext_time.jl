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

pInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_P.csv";
kInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_K.csv";
ϕInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_Phi_full.csv";

pData = readInP(pInputAdd,kInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
disData,Ω = autoUGen("LogNormal",[log(35),0.5],nameD,dparams,500,1 - pData.p0);
disData = orderdisData(disData,Ω);

text,xext,fext,gext,mp = extForm_cheat(pData,disData,Ω,999999);
ubmp = mp.objVal;
lbmp = mp.objBound;
gap = (mp.objVal - mp.objBound)/mp.objVal;
dDict = [text,xext,fext,gext,ubmp,lbmp,gap];
save("test_Ext_time.jld","dDict",dDict);
