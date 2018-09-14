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
@everywhere include("expModel.jl");
@everywhere include("ubCalFunc.jl");
@everywhere include("tighten.jl");
@everywhere include("partition_LP.jl");
@everywhere include("partition_LR.jl");
@everywhere include("part_tight.jl");
@everywhere include("partitionSolve.jl");
@everywhere include("convexify.jl");

pInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_P.csv";
kInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_K.csv";
ϕInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_Phi_full.csv";

pData = readInP(pInputAdd,kInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
disData,Ω = autoUGen("LogNormal",[log(35),0.5],nameD,dparams,5,1 - pData.p0);
disData = orderdisData(disData,Ω);

# deterministic solution
tdet,xdet,fdet = detBuild(pData);
ubdet = ubCalP(pData,disData,Ω,xdet,tdet,999999);

# expected solution
eH = exp(log(35)+1/2*0.5^2);
ed = Dict();
for i in pData.II
    if i == 0
        ed[i] = 0;
    else
        ed[i] = dparams[i][1];
    end
end
texp,xexp,fexp,Gexp,mexp = expModel(pData,eH,ed);
ubexp = ubCalP(pData,disData,Ω,xexp,texp,999999);

# our decomposition method
tic();
tbest,xbest,lbCost,ubCost = partitionSolve(pData,disData,0.001);
timedecomp = toc();
gapdecomp = (ubCost - lbCost)/ubCost;

# extensive formulation
tic();
text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,1e-2,999999);
timeext = toc();
θext = Dict();
for ω in Ω
    θext[ω] = getvalue(mp[:t][0,ω]);
end
ubmp = mext.objVal;
lbmp = mext.objBound;
gapext = (mext.objVal - mext.objBound)/mext.objVal;
ubext = ubCal(pData,disData,Ω,xext,text,999999);
dDict = [tdet,xdet,fdet,ubdet,texp,xexp,fexp,Gexp,ubexp,
            tbest,xbest,lbCost,ubCost,gapdecomp,timedecomp,
            text,xext,fext,gext,ubmp,lbmp,gapext,timeext];
save("test_Ext_time.jld","dDict",dDict);

############################################
ubList = [];
for n in 1:30
    disData,Ω = autoUGen("LogNormal",[log(35),0.5],nameD,dparams,500,1 - pData.p0);
    disData = orderdisData(disData,Ω);
    ubdet = ubCal(pData,disData,Ω,xdet,tdet,999999);
    ubexp = ubCal(pData,disData,Ω,xexp,texp,999999);
    ubext = ubCal(pData,disData,Ω,xext,text,999999);
    push!(ubList,[ubdet,ubexp,ubext]);
end
