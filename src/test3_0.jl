# test the variance of each model
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
nameD1,dparams1 = readInUnc(ϕInputAdd_1);
nameD2,dparams2 = readInUnc(ϕInputAdd_2);
nameD3,dparams3 = readInUnc(ϕInputAdd_3);
nameD4,dparams4 = readInUnc(ϕInputAdd_4);
println("------------- Finished Input Process -------------");

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
tList = [];
push!(tList,tdetO);
xList = [];
push!(xList,xdetO);

# sample N different points
N = 10;
ln = 1;
ad = Uniform(0,1);
while ln <= N
    xn = rand(ad,length(pData.II));
    if sum(xn) <= pData.B
        xnd = Dict();
        iin = 1;
        for i in pData.II
            xnd[i,1] = xn[iin];
            iin += 1;
        end
        push!(xList,xnd);
        tn = floyd(pData,xnd);
        push!(tList,tn);
        ln += 1;
    end
end

dDict = Dict();

disData1,Ω1 = autoUGen("Uniform",[35.125,35.125],nameD1,dparams1,5000);
σmax = -Inf;
for iSample in 1:11
    f1List = [];
    for ω in Ω1
        cω = subInt(pData,disData1[ω],xList[iSample],tList[iSample]);
        push!(f1List,cω);
    end
    if sum((item - mean(f1List)^2) for item in f1List)/length(f1List) > σmax
        σmax = sum((item - mean(f1List)^2) for item in f1List)/length(f1List);
    end
end
dDict[1] = σmax;

disData2,Ω2 = autoUGen("LogNormal",[log(35),0.5],nameD2,dparams2,5000);
σmax = -Inf;
for iSample in 1:11
    f2List = [];
    for ω in Ω2
        cω = subInt(pData,disData2[ω],xdet,tdet);
        push!(f2List,cω);
    end
    if sum((item - mean(f2List)^2) for item in f2List)/length(f2List) > σmax
        σmax = sum((item - mean(f2List)^2) for item in f2List)/length(f2List);
    end
end
dDict[2] = σmax;

f4List = [];
disData4,Ω4 = autoUGen("LogNormal",[log(35),0.5],nameD4,dparams4,5000);
σmax = -Inf;
for iSample in 1:11
    f4List = [];
    for ω in Ω2
        cω = subInt(pData,disData2[ω],xdet,tdet);
        push!(f4List,cω);
    end
    if sum((item - mean(f4List)^2) for item in f4List)/length(f4List) > σmax
        σmax = sum((item - mean(f4List)^2) for item in f4List)/length(f4List);
    end
end
dDict[4] = σmax;

save("test3.jld","dDict",dDict);
