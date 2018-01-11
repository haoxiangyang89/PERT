# serial example to test out the cut structure
@everywhere using JuMP,Gurobi,CPLEX,Cbc,Clp;

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

pInputAdd = "test_Serial_P.csv";
kInputAdd = "test_Serial_K.csv";
ΩInputAdd = "test_Serial_Omega.csv";

pData = readInP(pInputAdd,kInputAdd);
disData,Ω = readInDis(ΩInputAdd);

text,xext,fext,mext = extForm_cheat(pData,disData,Ω);
