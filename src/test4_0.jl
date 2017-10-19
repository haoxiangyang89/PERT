# test 4: test the solution of each model
@everywhere using JuMP,Gurobi,Cbc,Clp;

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

nameDF,dparamsF = readInUnc(ϕInputAdd_4);
disDataF,ΩF = autoUGen("LogNormal",[log(35),0.5],nameDF,dparamsF,10000);

pData = readInP(pInputAdd,kInputAdd);

# deterministic
tdet,xdet,fdet = detBuild(pData);
f0List = [];
for ω in ΩF
    cω = subInt(pData,disDataF[ω],xdet,tdet);
    push!(f0List,cω);
end
println("-------------- det Solved --------------")
# dOnly
nameD1,dparams1 = readInUnc(ϕInputAdd_1);
disData1,Ω1 = autoUGen("Uniform",[35.125,35.125],nameD1,dparams1,200);
text1,xext1,fext1,mext1 = extForm(pData,disData1,Ω1);
f1List = [];
for ω in ΩF
    cω = subInt(pData,disDataF[ω],xext1,text1);
    push!(f1List,cω);
end
println("-------------- dOnly Solved --------------");
# HOnly
nameD2,dparams2 = readInUnc(ϕInputAdd_2);
disData2,Ω2 = autoUGen("LogNormal",[log(35),0.5],nameD2,dparams2,200);
text2,xext2,fext2,mext2 = extForm(pData,disData2,Ω2);
f2List = [];
for ω in ΩF
    cω = subInt(pData,disDataF[ω],xext2,text2);
    push!(f2List,cω);
end
println("-------------- HOnly Solved --------------");
# Fixed
nameD3,dparams3 = readInUnc(ϕInputAdd_3);
disData3,Ω3 = autoUGen("Uniform",[35.125,35.125],nameD3,dparams3,1);
text3,xext3,fext3,mext3 = extForm(pData,disData3,Ω3);
f3List = [];
for ω in ΩF
    cω = subInt(pData,disDataF[ω],xext3,text3);
    push!(f3List,cω);
end
println("-------------- Fixed Solved --------------");
# Full
nameD4,dparams4 = readInUnc(ϕInputAdd_4);
disData4,Ω4 = autoUGen("LogNormal",[log(35),0.5],nameD4,dparams4,200);
text4,xext4,fext4,mext4 = extForm(pData,disData4,Ω4);
f4List = [];
for ω in ΩF
    cω = subInt(pData,disDataF[ω],xext4,text4);
    push!(f4List,cω);
end
println("-------------- Full Solved --------------");
save("test4_f1.jld","f0List",f0List,"f1List",f1List,"f2List",f2List,"f3List",f3List,"f4List",f4List);
save("test4_xt1.jld","x0",xdet,"x1",xext1,"x2",xext2,"x3",xext3,"x4",xext4,"t0",tdet,"t1",text1,"t2",text2,"t3",text3,"t4",text4);
