# test formulation where the disruption can impact the activity after its start
using Distributed;
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,Ipopt,MathProgBase;
@everywhere using Distributions,HDF5,JLD,DelimitedFiles,Statistics,SharedArrays;
@everywhere const GUROBI_ENV = Gurobi.Env();

# test sbb
@everywhere include("header.jl");
@everywhere include("pullDecomp.jl");

pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
            "/home/haoxiang/scratch/PERT_tests/current/14/",
            "/home/haoxiang/scratch/PERT_tests/current/19/",
            "/home/haoxiang/scratch/PERT_tests/current/35/",
            "/home/haoxiang/scratch/PERT_tests/current/55/",
            "/home/haoxiang/scratch/PERT_tests/current/75/"];
# pathList = ["/scratch/haoxiang/current/11/",
#                "/scratch/haoxiang/current/14/",
#                "/scratch/haoxiang/current/19/",
#                "/scratch/haoxiang/current/35/",
#                "/scratch/haoxiang/current/55/",
#                "/scratch/haoxiang/current/75/"];

Ωsize = [100,200,500,1000];
dDict = Dict();
ϵ = 1e-2;
sNList = [10,20,25,25];
MMList = [10,10,20,40];

for fileInd in 1:4
   filePath = pathList[fileInd];
   # compile the functions
   dDict[fileInd] = Dict();
   for Ωl in 1:length(Ωsize)
      global Ω = 1:Ωsize[Ωl];
      randNo = 1;
      sN = sNList[Ωl];
      MM = MMList[Ωl];
      extBool = true;

      # test the after-start disruption effect
      pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1,1,"test_P.csv","test_K.csv","test_Phi.csv","test_H.csv",
                      0, 0, 0, true);
      global pData = pData;
      disDataRaw = load(pathList[fileInd]*"solData_$(Ωsize[Ωl]).jld");
      disData = disDataRaw["data"][randNo];

      global allSucc = findSuccAll(pData);
      global distanceDict = Dict();
      for i in pData.II
         for j in allSucc[i]
             distanceDict[i,j] = detCal(pData,i,j);
         end
      end
      # text_a,xext_a,fext_a,gext_a,mext_a,timeext_a = extForm_cheat_after(pData,disData,Ω,ϵ,999999);
      # ubext_a = mext_a.objVal;
      # lbext_a = mext_a.objBound;
      # gapext_a = (ubext_a - lbext_a)/ubext_a;
      tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_after(pData,disData,Ω,sN,MM,noThreads,5,6,5,1e-2,5,21600);
      gapdecomp = (ubFull - lbFull)/ubFull;
      dDict[fileInd][Ωsize[Ωl]] = [tFull,xFull,ubFull,lbFull,gapdecomp,timedecomp];
      save("test_Ext_after.jld","dDict",dDict);
   end
end
