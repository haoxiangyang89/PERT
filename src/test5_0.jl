# test the discrete distribution with everything up to change, compare the results
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
ϕInputAdd_4 = "test_14_Phi_full.csv";
pData = readInP(pInputAdd,kInputAdd);
nameDF,dparamsF = readInUnc(ϕInputAdd_4);

disDataF,ΩF = autoUGen("LogNormal",[log(35),0.5],nameDF,dparamsF,200);

tdet,xdet,fdet = detBuild(pData);
tdetO = Dict();
xdetO = Dict();
for i in pData.II
    tdetO[i] = tdet[i];
    for j in pData.Ji[i]
        xdetO[i,j] = xdet[i,j];
    end
end
fdetU = ubCal(pData,disDataF,ΩF,xdet,tdet);
println("-------------- det Solved --------------");

tfull,xfull,ffull,mfull = extForm(pData,disDataF,ΩF);
tfullO = Dict();
xfullO = Dict();
for i in pData.II
    tfullO[i] = tfull[i];
    for j in pData.Ji[i]
        xfullO[i,j] = xfull[i,j];
    end
end
println("-------------- Full Solved --------------");

meanH = 0;
meand = Dict();
for i in pData.II
    if i != 0
        meand[i] = 0;
    end
end
for ω in ΩF
    meanH += disDataF[ω].H;
    for i in pData.II
        if i != 0
            meand[i] += disDataF[ω].d[i];
        end
    end
end
meanH = meanH/length(ΩF);
for i in pData.II
    if i != 0
        meand[i] = meand[i]/length(ΩF);
    end
end


disData1 = copy(disDataF);
for ω in ΩF
    disData1[ω].H = meanH;
end
tdonly,xdonly,fdonly,mdonly = extForm(pData,disData1,ΩF);
tdonlyO = Dict();
xdonlyO = Dict();
for i in pData.II
    tdonlyO[i] = tdonly[i];
    for j in pData.Ji[i]
        xdonlyO[i,j] = xdonly[i,j];
    end
end
fdonlyU = ubCal(pData,disDataF,ΩF,xdonly,tdonly);
println("-------------- dOnly Solved --------------");

disData2 = copy(disDataF);
for ω in ΩF
    disData2[ω].d = meand;
end
tHonly,xHonly,fHonly,mHonly = extForm(pData,disData2,ΩF);
tHonlyO = Dict();
xHonlyO = Dict();
for i in pData.II
    tHonlyO[i] = tHonly[i];
    for j in pData.Ji[i]
        xHonlyO[i,j] = xHonly[i,j];
    end
end
fHonlyU = ubCal(pData,disDataF,ΩF,xHonly,tHonly);
println("-------------- HOnly Solved --------------");

disData3 = copy(disData2);
for ω in ΩF
    disData3[ω].H = meanH;
end
tfixed,xfixed,ffixed,mfixed = extForm(pData,disData3,ΩF);
tfixedO = Dict();
xfixedO = Dict();
for i in pData.II
    tfixedO[i] = tfixed[i];
    for j in pData.Ji[i]
        xfixedO[i,j] = xfixed[i,j];
    end
end
ffixedU = ubCal(pData,disDataF,ΩF,xfixed,tfixed);
println("-------------- Fixed Solved --------------");
save("test5_f1.jld","fdet",fdetU,"ffull",ffull,"fdonly",fdonlyU,"fHonly",fHonlyU,"ffixed",ffixedU);
save("test5_xt1.jld","xdet",xdetO,"xfull",xfullO,"xdonly",xdonlyO,"xHonly",xHonlyO,"xfixed",xfixedO,
    "tdet",tdetO,"tfull",tfullO,"tdonly",tdonlyO,"tHonly",tHonlyO,"tfixed",tfixedO);
