# This is the script to load the data and solve the problem
@everywhere using JuMP,Gurobi;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("main_1l.jl");
@everywhere include("cutGen_1l.jl");
@everywhere include("tightenGen_1l.jl");
@everywhere include("branchFunc.jl");
@everywhere include("detFunc_1l.jl");
@everywhere include("extForm_1l.jl");

pInputAdd = "test_Full_P_3.csv";
kInputAdd = "test_Full_K.csv";
ΩInputAdd = "test_Full_Omega_3.csv";

pData = readInP(pInputAdd,kInputAdd);
disData,Ω = readInDis(ΩInputAdd);

@time tbest,xbest,fbest = cutProc_Benders(pData,disData,Ω);
@time text,xext,fext,mext = extForm(pData,disData,Ω);
