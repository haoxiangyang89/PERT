using Distributed;
addprocs(20);
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

#filePath = "/Users/haoxiangyang/Desktop/PERT_tests/75_Lognormal_Exponential/"
filePath = "/home/haoxiang/PERT_tests/14_Lognormal_Exponential";
Ωsize = 500;
Ω = 1:Ωsize;
ϵ = 1e-2;
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize);
data = load("test_cuts.jld");
disData = data["disData"];
dataDet = load("test_Ext_time_exponential.jld");
allSucc = findSuccAll(pData);
distanceDict = Dict();
for i in pData.II
    for j in allSucc[i]
        distanceDict[i,j] = detCal(pData,i,j);
    end
end
#disData = disDataSet[1];

# deterministic solution
tdet,xdet,fdet = detBuild(pData);
ubdet = ubCalP(pData,disData,Ω,xdet,tdet,999999);

# expected solution
eH = mean(buildDistrn(nameH,Hparams));
ed = Dict();
for i in pData.II
    if i == 0
        ed[i] = 0;
    else
        ed[i] = mean(buildDistrn(nameD,dparams[i]));
    end
end
texp,xexp,fexp,Gexp,mexp = expModel(pData,eH,ed);
ubexp = ubCalP(pData,disData,Ω,xexp,texp,999999);

# our decomposition method
global sN = 20;
global MM = 25;
tempTimer = time();
include("partSolve_Callback_tightened_sol.jl");
timedecomp = time() - tempTimer;
gapdecomp = (ubCost - lbCost)/ubCost;
ubFull = ubCost;
lbFull = lbCost;
xFull = deepcopy(xbest);
tFull = deepcopy(tbest);

# extensive formulation
tempTimer = time();
text,xext,fext,gext,mext = extForm_cheat(pData,disData,Ω,1e-4,999999);
timeext = time() - tempTimer;
θext = Dict();
for ω in Ω
    θext[ω] = getvalue(mext[:t][0,ω]);
end
ubmp = mext.objVal;
lbmp = mext.objBound;
gapext = (mext.objVal - mext.objBound)/mext.objVal;
ubext = ubCal(pData,disData,Ω,xext,text,999999);
print(ubdet," ",ubexp," ",ubext);
dDict = [tdet,xdet,fdet,ubdet,texp,xexp,fexp,Gexp,ubexp,
            tFull,xFull,lbFull,ubFull,gapdecomp,timedecomp];
            text,xext,fext,gext,ubmp,lbmp,gapext,timeext];
save("test_Ext_time_d.jld","dDict",dDict);

############################################
ubList = [];
for n in 1:30
    disData1,Ω = autoUGen("LogNormal",Hparams,"Exponential",dparams,500,1 - pData.p0);
    disData1 = orderdisData(disData1,Ω);
    ubdet1 = ubCalP(pData,disData1,Ω,xdet,tdet,999999);
    ubexp1 = ubCalP(pData,disData1,Ω,xexp,texp,999999);
    ubbest1 = ubCalP(pData,disData1,Ω,xbest,tbest,999999);
    push!(ubList,[ubdet1,ubexp1,ubbest1]);
end
