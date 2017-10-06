# This is the script to load the data and solve the problem
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
ΩInputAdd = "test_Full_Omega.csv";
ϕInputAdd = "test_Full_Phi.csv";

pData = readInP(pInputAdd,kInputAdd);
#disData,Ω = readInDis(ΩInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
disData,Ω = autoUGen("LogNormal",[log(35),0.5],nameD,dparams,100);

@time tbest,xbest,fbest = cutProc_Benders(pData,disData,Ω);
@time text,xext,fext,mext = extForm(pData,disData,Ω);

# small test for consistency
# ezn = 0;
# for i in 1:10
#     disData,Ω = autoUGen("LogNormal",[log(35),0.5],nameD,dparams,10);
#     @time text,xext,fext,mext = extForm(pData,disData,Ω);
#     ezn += fext;
# end
# ezn = ezn/10;
