@everywhere using JuMP,Gurobi,CPLEX,Cbc,Clp,HDF5,JLD;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("main_1l.jl");
@everywhere include("main_pull.jl");
@everywhere include("cutGen_1l.jl");
@everywhere include("tightenGen_1l.jl");
@everywhere include("branchFunc.jl");
@everywhere include("detFunc_1l.jl");
@everywhere include("extForm_1l.jl");
@everywhere include("main_pull.jl");

pInputAdd = "test_14_P.csv";
kInputAdd = "test_14_K.csv";
ΩInputAdd = "test_14_Omega_full.csv";
ϕInputAdd = "test_14_Phi_full.csv";

pData = readInP(pInputAdd,kInputAdd);
#disData,Ω = readInDis(ΩInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);

vList1 = [];
vList2 = [];
vList3 = [];
vList4 = [];
vList5 = [];
vList6 = [];

for i in 1:2000
    disData,Ω = autoUGenStrata2("LogNormal",[log(35),0.5],nameD,dparams,1,50,1-pData.p0);
    ωSeq,ωDict,Hω,Ωt = getOmegaSeq(disData);
    text,xext,fext,mext = extForm(pData,disData,Ω);
    push!(vList1,fext);

    disData,Ω = autoUGenStrata2("LogNormal",[log(35),0.5],nameD,dparams,2,25,1-pData.p0);
    ωSeq,ωDict,Hω,Ωt = getOmegaSeq(disData);
    text,xext,fext,mext = extForm(pData,disData,Ω);
    push!(vList2,fext);

    disData,Ω = autoUGenStrata2("LogNormal",[log(35),0.5],nameD,dparams,5,10,1-pData.p0);
    ωSeq,ωDict,Hω,Ωt = getOmegaSeq(disData);
    text,xext,fext,mext = extForm(pData,disData,Ω);
    push!(vList3,fext);

    disData,Ω = autoUGenStrata2("LogNormal",[log(35),0.5],nameD,dparams,10,5,1-pData.p0);
    ωSeq,ωDict,Hω,Ωt = getOmegaSeq(disData);
    text,xext,fext,mext = extForm(pData,disData,Ω);
    push!(vList4,fext);

    disData,Ω = autoUGenStrata2("LogNormal",[log(35),0.5],nameD,dparams,25,2,1-pData.p0);
    ωSeq,ωDict,Hω,Ωt = getOmegaSeq(disData);
    text,xext,fext,mext = extForm(pData,disData,Ω);
    push!(vList5,fext);

    disData,Ω = autoUGenStrata2("LogNormal",[log(35),0.5],nameD,dparams,50,1,1-pData.p0);
    ωSeq,ωDict,Hω,Ωt = getOmegaSeq(disData);
    text,xext,fext,mext = extForm(pData,disData,Ω);
    push!(vList6,fext);
    println("------------------ Iteration $(i) Finished ------------------");
end

save("consistency_2.jld","v1",vList1,"v2",vList2,"v3",vList3,"v4",vList4,"v5",vList5,"v6",vList6);
