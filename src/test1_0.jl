# This is the script to load the data and solve the problem
@everywhere using JuMP,Gurobi;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("main_1l.jl");
@everywhere include("cutGen_1l.jl");
@everywhere include("tightenGen_1l.jl");
@everywhere include("branchFunc.jl");

pInputAdd = "test_full_P.csv";
kInputAdd = "test_full_K.csv";
ΩInputAdd = "test_full_Omega.csv";

pData = readInP(pInputAdd,kInputAdd);
disData,Ω = readInDis(ΩInputAdd);

tbest,xbest,fbest = cutProc_Benders(pData,disData,Ω);
