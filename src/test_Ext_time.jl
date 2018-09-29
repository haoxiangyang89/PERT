addprocs(20);
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

#filePath = "/Users/haoxiangyang/Desktop/PERT_tests/14_ExponentialD_LogNormalH/"
filePath = "/home/haoxiang/PERT_tests/14";
Ωsize = 500;
Ω = 1:Ωsize;
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize);
disData = disDataSet[1];

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
tbest,xbest,lbCost,ubCost = partitionSolve(pData,disData,0.01);
timedecomp = toc();
gapdecomp = (ubCost - lbCost)/ubCost;

# extensive formulation
tic();
text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,1e-4,999999);
timeext = toc();
θext = Dict();
for ω in Ω
    θext[ω] = getvalue(mext[:t][0,ω]);
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
    disData1,Ω = autoUGen("LogNormal",Hparams,"Exponential",dparams,500,1 - pData.p0);
    disData1 = orderdisData(disData1,Ω);
    ubdet1 = ubCalP(pData,disData1,Ω,xdet,tdet,999999);
    ubexp1 = ubCalP(pData,disData1,Ω,xexp,texp,999999);
    ubext1 = ubCalP(pData,disData1,Ω,xext,text,999999);
    push!(ubList,[ubdet1,ubexp1,ubext1]);
end
