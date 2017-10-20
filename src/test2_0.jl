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
nameD1,dparams1 = readInUnc(ϕInputAdd_1);
nameD2,dparams2 = readInUnc(ϕInputAdd_2);
nameD3,dparams3 = readInUnc(ϕInputAdd_3);
nameD4,dparams4 = readInUnc(ϕInputAdd_4);
println("------------- Finished Input Process -------------");

dDict = Dict();

# only magnitude changes
text1List = [];
xext1List = [];
fext1List = [];
for l in 1:10
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
    push!(text1List,text1O);
    push!(xext1List,xext1O);
    push!(fext1List,fext1);
    println("------------- dOnly No. $(l) finished -------------");
end
dDict[1] = (text1List,xext1List,fext1List);

# only time changes
text2List = [];
xext2List = [];
fext2List = [];
for l in 1:10
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
    push!(text2List,text2O);
    push!(xext2List,xext2O);
    push!(fext2List,fext2);
    println("------------- HOnly No. $(l) finished -------------");
end
dDict[2] = (text1List,xext1List,fext1List);

# fixed
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
dDict[3] = (text3O,xext3O,fext3);

# both
text4List = [];
xext4List = [];
fext4List = [];
for l in 1:10
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
    push!(text4List,text4O);
    push!(xext4List,xext4O);
    push!(fext4List,fext4);
    println("------------- Full No. $(l) finished -------------");
end
dDict[4] = (text4List,xext4List,fext4List);

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
dDict[5] = (tdetO,xdetO,fdet);
save("test2.jld","dDict",dDict);
