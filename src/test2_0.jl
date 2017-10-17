# test the variance of the recourse function with the deterministic solution
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

# only magnitude changes
nameD1,dparams1 = readInUnc(ϕInputAdd_1);
text1List = [];
xext1List = [];
fext1List = [];
for i in 1:10
    disData1,Ω1 = autoUGen("Uniform",[35.125,35.125],nameD1,dparams1,100);
    text1,xext1,fext1,mext1 = extForm(pData,disData1,Ω1);
    push!(text1List,text1);
    push!(xext1List,xext1);
    push!(xext1List,fext1);
end

# only time changes
nameD2,dparams2 = readInUnc(ϕInputAdd_2);
text2List = [];
xext2List = [];
fext2List = [];
for i in 1:10
    disData2,Ω2 = autoUGen("LogNormal",[log(35),0.5],nameD2,dparams2,100);
    text2,xext2,fext2,mext2 = extForm(pData,disData2,Ω2);
    push!(text2List,text2);
    push!(xext2List,xext2);
    push!(xext2List,fext2);
end

# fixed
nameD3,dparams3 = readInUnc(ϕInputAdd_3);
disData3,Ω3 = autoUGen("Uniform",[35.125,35.125],nameD3,dparams3,1);
text3,xext3,fext3,mext3 = extForm(pData,disData3,Ω3);

# both
nameD4,dparams4 = readInUnc(ϕInputAdd_4);
text4List = [];
xext4List = [];
fext4List = [];
for i in 1:10
    disData4,Ω4 = autoUGen("LogNormal",[log(35),0.5],nameD4,dparams4,100);
    text4,xext4,fext4,mext4 = extForm(pData,disData4,Ω4);
    push!(text4List,text4);
    push!(xext4List,xext4);
    push!(xext4List,fext4);
end

# deterministic
tdet,xdet,fdet = detBuild(pData);
