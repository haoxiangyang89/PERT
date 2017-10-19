# test the variance of the recourse function with the deterministic solution
@everywhere using JuMP,Gurobi,Cbc,Clp,HDF5,JLD;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("main_1l.jl");
@everywhere include("cutGen_1l.jl");
@everywhere include("tightenGen_1l.jl");
@everywhere include("branchFunc.jl");
@everywhere include("detFunc_1l.jl");
@everywhere include("extForm_1l.jl");

pInputAdd = "test_14_P.csv";
kInputAdd = "test_14_K.csv";
ϕInputAdd_1 = "test_14_Phi_dOnly.csv";
ϕInputAdd_2 = "test_14_Phi_HOnly.csv";
ϕInputAdd_3 = "test_14_Phi_fixed.csv";
ϕInputAdd_4 = "test_14_Phi_full.csv";

pData = readInP(pInputAdd,kInputAdd);

dDict = Dict();

# only magnitude changes
nameD1,dparams1 = readInUnc(ϕInputAdd_1);
text1List = [];
xext1List = [];
fext1List = [];
for i in 1:10
    disData1,Ω1 = autoUGen("Uniform",[35.125,35.125],nameD1,dparams1,200);
    text1,xext1,fext1,mext1 = extForm(pData,disData1,Ω1);
    push!(text1List,text1);
    push!(xext1List,xext1);
    push!(fext1List,fext1);
    println("------------- dOnly No. $(i) finished -------------");
end
dDict[1] = (text1List,xext1List,fext1List);

# only time changes
nameD2,dparams2 = readInUnc(ϕInputAdd_2);
text2List = [];
xext2List = [];
fext2List = [];
for i in 1:10
    disData2,Ω2 = autoUGen("LogNormal",[log(35),0.5],nameD2,dparams2,200);
    text2,xext2,fext2,mext2 = extForm(pData,disData2,Ω2);
    push!(text2List,text2);
    push!(xext2List,xext2);
    push!(fext2List,fext2);
    println("------------- HOnly No. $(i) finished -------------");
end
dDict[2] = (text1List,xext1List,fext1List);

# fixed
nameD3,dparams3 = readInUnc(ϕInputAdd_3);
disData3,Ω3 = autoUGen("Uniform",[35.125,35.125],nameD3,dparams3,1);
text3,xext3,fext3,mext3 = extForm(pData,disData3,Ω3);
dDict[3] = (text3,xext3,fext3);

# both
nameD4,dparams4 = readInUnc(ϕInputAdd_4);
text4List = [];
xext4List = [];
fext4List = [];
for i in 1:10
    disData4,Ω4 = autoUGen("LogNormal",[log(35),0.5],nameD4,dparams4,200);
    text4,xext4,fext4,mext4 = extForm(pData,disData4,Ω4);
    push!(text4List,text4);
    push!(xext4List,xext4);
    push!(fext4List,fext4);
    println("------------- Full No. $(i) finished -------------");
end
dDict[4] = (text4List,xext4List,fext4List);

# deterministic
tdet,xdet,fdet = detBuild(pData);
dDict[5] = (tdet,xdet,fdet);
save("test2.jld","dDict",dDict);
