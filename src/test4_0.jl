# test 4: test the solution of each model
using Distributed;
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

pData = readInP(pInputAdd,kInputAdd);
nameDF,dparamsF = readInUnc(ϕInputAdd_4);
nameD1,dparams1 = readInUnc(ϕInputAdd_1);
nameD2,dparams2 = readInUnc(ϕInputAdd_2);
nameD3,dparams3 = readInUnc(ϕInputAdd_3);
nameD4,dparams4 = readInUnc(ϕInputAdd_4);
println("------------- Finished Input Process -------------");

disDataF,ΩF = autoUGen("LogNormal",[log(35),0.5],nameDF,dparamsF,10000);

# deterministic
tdet,xdet,fdet = detBuild(pData);
tdetO = Dict();
xdetO = Dict();
for i in pData.II
    tdetO[i] = tdet[i];
    for j in pData.Ji[i]
        xdetO[i,j] = xdet[i,j];
    end
end
f0List = [];
for ω in ΩF
    cω = subInt(pData,disDataF[ω],xdet,tdet);
    push!(f0List,cω);
end
println("-------------- det Solved --------------")
# dOnly
disData1,Ω1 = autoUGen("Uniform",[35.125,35.125],nameD1,dparams1,200);
text1,xext1,fext1,mext1 = extForm(pData,disData1,Ω1);
text1O = Dict();
xext1O = Dict();
for i in pData.II
    text1O[i] = text1[i];
    for j in pData.Ji[i]
        xext1O[i,j] = xext1[i,j];
    end
end
f1List = [];
for ω in ΩF
    cω = subInt(pData,disDataF[ω],xext1,text1);
    push!(f1List,cω);
end
println("-------------- dOnly Solved --------------");
# HOnly
disData2,Ω2 = autoUGen("LogNormal",[log(35),0.5],nameD2,dparams2,200);
text2,xext2,fext2,mext2 = extForm(pData,disData2,Ω2);
text2O = Dict();
xext2O = Dict();
for i in pData.II
    text2O[i] = text2[i];
    for j in pData.Ji[i]
        xext2O[i,j] = xext2[i,j];
    end
end
f2List = [];
for ω in ΩF
    cω = subInt(pData,disDataF[ω],xext2,text2);
    push!(f2List,cω);
end
println("-------------- HOnly Solved --------------");
# Fixed
disData3,Ω3 = autoUGen("Uniform",[35.125,35.125],nameD3,dparams3,1);
text3,xext3,fext3,mext3 = extForm(pData,disData3,Ω3);
text3O = Dict();
xext3O = Dict();
for i in pData.II
    text3O[i] = text3[i];
    for j in pData.Ji[i]
        xext3O[i,j] = xext3[i,j];
    end
end
f3List = [];
for ω in ΩF
    cω = subInt(pData,disDataF[ω],xext3,text3);
    push!(f3List,cω);
end
println("-------------- Fixed Solved --------------");
# Full
disData4,Ω4 = autoUGen("LogNormal",[log(35),0.5],nameD4,dparams4,200);
text4,xext4,fext4,mext4 = extForm(pData,disData4,Ω4);
text4O = Dict();
xext4O = Dict();
for i in pData.II
    text4O[i] = text4[i];
    for j in pData.Ji[i]
        xext4O[i,j] = xext4[i,j];
    end
end
f4List = [];
for ω in ΩF
    cω = subInt(pData,disDataF[ω],xext4,text4);
    push!(f4List,cω);
end
println("-------------- Full Solved --------------");
save("test4_f1.jld","f0List",f0List,"f1List",f1List,"f2List",f2List,"f3List",f3List,"f4List",f4List);
save("test4_xt1.jld","x0",xdetO,"x1",xext1O,"x2",xext2O,"x3",xext3O,"x4",xext4O,"t0",tdetO,"t1",text1O,"t2",text2O,"t3",text3O,"t4",text4O);
