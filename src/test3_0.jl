# test the variance of each model

@everywhere using JuMP,Gurobi,Cbc,Clp;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("main_1l.jl");
@everywhere include("cutGen_1l.jl");
@everywhere include("tightenGen_1l.jl");
@everywhere include("branchFunc.jl");
@everywhere include("detFunc_1l.jl");
@everywhere include("extForm_1l.jl");

pInputAdd = "test_Full_P.csv";
kInputAdd = "test_Full_K.csv";
ϕInputAdd_1 = "test_Full_Phi_dOnly.csv";
ϕInputAdd_2 = "test_Full_Phi_HOnly.csv";
ϕInputAdd_3 = "test_Full_Phi_fixed.csv";
ϕInputAdd_4 = "test_Full_Phi_full.csv";

pData = readInP(pInputAdd,kInputAdd);

# deterministic
tdet,xdet,fdet = detBuild(pData);

nameD1,dparams1 = readInUnc(ϕInputAdd_1);
f1List = [];
disData1,Ω1 = autoUGen("Uniform",[35.125,35.125],nameD1,dparams1,1000);
for ω in Ω1
    cω = subInt(pData,disData1[ω],xdet,tdet);
    push!(f1List,cω);
end

nameD2,dparams2 = readInUnc(ϕInputAdd_2);
f2List = [];
disData2,Ω2 = autoUGen("LogNormal",[log(35),0.5],nameD2,dparams2,1000);
for ω in Ω2
    cω = subInt(pData,disData2[ω],xdet,tdet);
    push!(f2List,cω);
end

nameD4,dparams4 = readInUnc(ϕInputAdd_4);
f4List = [];
disData4,Ω4 = autoUGen("LogNormal",[log(35),0.5],nameD4,dparams4,1000);
for ω in Ω4
    cω = subInt(pData,disData4[ω],xdet,tdet);
    push!(f4List,cω);
end
