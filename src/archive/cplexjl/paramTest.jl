# parametrized test
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

##==============================================================================
# set up the test here, changing the parameters (serial network type 1)
n = 2;
k = 1000;
D0 = 1;
d0 = 1.99;
e0 = 0.9;
p0 = 0.5;
B = 1;
b0 = 1;
nameH = "LogNormal";
Hparams = [2,0.5];
Ωsize = 10;

pData,disData,Ω = serialMaker1(n,k,D0,d0,e0,p0,B,b0,nameH,Hparams,Ωsize);

distrD = Dict();
distrName = Dict();
for i in 1:n
    distrName[i] = "Singleton";
    distrD[i] = disData[1].d[i]
end

##==============================================================================
# generic network
n = 14;
D0 = [5,30,25,20,15,24,30,49,40,30,45,25,21,5];
B = 3;
e0 = 1/2;
p0 = 0.99;
pd = 0.01;

pData,disData,Ω = serialMaker2(n,k,D0,d0,e0,p0,B,b0,nameH,Hparams,Ωsize);

##==============================================================================
# special case network
#n = 10;                 # how many parallel arcs in the bottom path
D0 = 1;
d0 = 1.99;
k = 10000;
k1 = 12000;
b0 = 1;
B = 2;
e0 = 0.9;
p0 = 0.5;
pd = 0.5;
nameH = "Uniform";
Hparams = [1.5,2.5];
Ωsize = 2;
distrD = Dict();
distrName = Dict();
for i in 0:5
    distrName[i] = "Singleton";
end
distrD[1] = d0;
distrD[2] = k*d0;
distrD[3] = 0;
distrD[4] = 0;
distrD[0] = 0;
distrD[5] = (k1-1)*D0 - pd;
distrName[6] = "Categorical";
distrD[6] = [[1-pd,pd],[0,(k1-1)*D0/pd - 1]];

pData,disData,Ω = specialMaker(k,k1,D0,d0,e0,p0,pd,B,b0,nameH,Hparams,Ωsize,distrName,distrD);

##==============================================================================
# generic test cases builder
II = [0,1,2,3,4];
k1 = 7000;
b0 = 1;
D0 = 1;
e0 = 0.99;
B = 1;
D = Dict(0 => 0,1 => k1*D0,2 => D0,3 => D0,4 => D0);
K = [(1,2),(2,3),(2,4),(3,0),(4,0)];
Pre = Dict();
Succ = Dict();
for i in II
    Pre[i] = [];
    Succ[i] = [];
end
Pre[2] = [1];
Pre[3] = [2];
Pre[4] = [2];
Succ[1] = [2];
Succ[2] = [3,4];

for i in II
    if i != 0
        push!(Pre[0],i);
        push!(Succ[i],0);
    end
end

p0 = 0.1;
pd = 0.5;
Ji = Dict();
Ji[0] = [];

b = Dict();
b[0] = Dict();
eff = Dict();
eff[0] = Dict();
for i in II
    Ji[i] = [1];
    b[i] = Dict();
    b[i][1] = b0;
    eff[i] = Dict();
    eff[i][1] = e0;
end

pData = pInfo(II,Ji,D,b,eff,B,p0,K,Pre,Succ);

nameH = "LogNormal";
Hparams = [2,0.5];
Ωsize = 10;
distrD = Dict();
distrName = Dict();
for i in 0:3
    distrName[i] = "Singleton";
    distrD[i] = 0;
end
distrD[3] = (k1-1)*D0 - pd;
distrName[4] = "Categorical";
distrD[4] = [[1-pd,pd],[0,(k1-1)*1/pd - 1]];
disData = Dict();
disData,Ω = autoUGen(nameH,Hparams,distrName,distrD,Ωsize,1 - pData.p0);
disData = orderdisData(disData,Ω);

##==============================================================================
# running tests codes
global ϵ = 1e-4;
pdData = deepcopy(pData);
for i in pData.II
    if i != 0
        pdData.D[i] = pData.D[i] + maximum([disData[ω].d[i] for ω in Ω])
    else
        pdData.D[i] = pData.D[i];
    end
end
lDict = longestPath(pdData);
for i in pData.II
    lDict[i] += disData[length(Ω)].H;
end
Tmax1 = lDict[0];

# deterministic solution
tdet,xdet,fdet = detBuild(pData);
ubdet = ubCalP(pData,disData,Ω,xdet,tdet,Tmax1);

# expected solution
eH = mean(buildDistrn(nameH,Hparams));
ed = Dict();
for i in pData.II
    if i == 0
        ed[i] = 0;
    else
        ed[i] = mean(buildDistrn(distrName[i],distrD[i]));
    end
end
texp,xexp,fexp,Gexp,mexp = expModel(pData,eH,ed,Tmax1);
ubexp = ubCalP(pData,disData,Ω,xexp,texp,Tmax1);

# dOnly solution
global disData1;
disData1 = deepcopy(disData);
for ω in Ω
    disData[ω].H = mean(buildDistrn(nameH,Hparams));
end
tdOnly,xdOnly,fdOnly,gdOnly,mdOnly = extForm_cheat(pData,disData,Ω,1e-4,Tmax1);
disData = deepcopy(disData1);
ubdOnly = ubCalP(pData,disData,Ω,xdOnly,tdOnly,Tmax1);

# HOnly solution
for ω in Ω
    for i in pData.II
        if i != 0
            disData[ω].d[i] = mean(buildDistrn(distrName[i],distrD[i]));
            if disData[ω].d[i] < 1e-4
                disData[ω].d[i] = 0;
            end
        end
    end
end
tHOnly,xHOnly,fHOnly,gHOnly,mHOnly = extForm_cheat(pData,disData,Ω,1e-4,Tmax1);
disData = deepcopy(disData1);
ubHOnly = ubCalP(pData,disData,Ω,xHOnly,tHOnly,Tmax1);

# full solution
#include("partSolve_Callback_tightened.jl");
# gapFull = (ubCost - lbCost)/ubCost;
# ubFull = ubCost;
# lbFull = lbCost;
# xFull = deepcopy(xbest);
# tFull = deepcopy(tbest);
tFull,xFull,fFull,gFull,mFull = extForm_cheat(pData,disData,Ω,1e-4,Tmax1);
ubFull = ubCalP(pData,disData,Ω,xFull,tFull,Tmax1);


dDict = [tdet,xdet,fdet,ubdet,
            texp,xexp,fexp,Gexp,ubexp,
            tFull,xFull,ubFull,lbFull,
            tdOnly,xdOnly,gdOnly,ubdOnly,
            tHOnly,xHOnly,ubHOnly];
