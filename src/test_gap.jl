# test the gap between the lower bound and the upper bound
addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
@everywhere const GUROBI_ENV = Gurobi.Env()

# test sbb
@everywhere include("header.jl");

pathList = ["/scratch/haoxiang/current/11/",
            "/scratch/haoxiang/current/14/",
            "/scratch/haoxiang/current/19/",
            "/scratch/haoxiang/current/35/",
            "/scratch/haoxiang/current/55/",
            "/scratch/haoxiang/current/75/"];
# pathList = ["/home/haoxiang/scratch/PERT_tests/current/11/",
#             "/home/haoxiang/scratch/PERT_tests/current/14/",
#             "/home/haoxiang/scratch/PERT_tests/current/19/",
#             "/home/haoxiang/scratch/PERT_tests/current/35/",
#             "/home/haoxiang/scratch/PERT_tests/current/55/",
#             "/home/haoxiang/scratch/PERT_tests/current/75/"];

dDict = Dict();
ubDict = Dict();
Ωsize = [10,20,50,100,200,500,1000];
sNList = [10,20,50,100,20,25,25];
MMList = [1,1,1,1,10,20,40];
global ϵ = 1e-2;

for fileInd in 1:4
    # load the data
    filePath = pathList[fileInd]
    simuPath = filePath*"simuData.jld";
    data5000Raw = load(simuPath);

    dDict[fileInd] = Dict();
    ubDict[fileInd] = Dict();
    pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1);
    global pData = pData;
    global allSucc = findSuccAll(pData);
    global distanceDict = Dict();
    for i in pData.II
        for j in allSucc[i]
            distanceDict[i,j] = detCal(pData,i,j);
        end
    end

    for Ωl in 1:length(Ωsize)
        global Ω = 1:Ωsize[Ωl];
        global sN = sNList[Ωl];
        global MM = MMList[Ωl];
        disData1 = data5000Raw["data"][1];
        Ω1 = 1:length(disData1);
        # generate a candidate solution by SAA
        ub5000List = [];
        tCanBest = Dict();
        xCanBest = Dict();
        for ii in 1:20
            pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize[Ωl]);
            disData = disDataSet[1];
            tCan,xCan,ubCan,lbCan,timeIterCan,treeListCan,timedecompCan = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,ϵ,5,10800);
            ub5000 = ubCalP(pData,disData1,Ω1,xCan,tCan,999999);
            push!(ub5000List,ub5000);
            if ub5000 <= minimum(ub5000List)
                tCanBest = deepcopy(tCan);
                xCanBest = deepcopy(xCan);
            end
        end
        ubDict[fileInd][Ωsize[Ωl]] = [tCanBest,xCanBest,ub5000List];

        #disDataRaw = load(pathList[fileInd]*"solData_$(Ωsize[Ωl]).jld");
        #disDataSet = disDataRaw["data"];
        pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize[Ωl],20000/Ωsize[Ωl]);
        dDict[fileInd][Ωsize[Ωl]] = [];
        n = 1;
        while n <= 20000/Ωsize[Ωl]
            # try
            global disData = disDataSet[n];

            # our decomposition method
            if Ωl <= 4
                tFull,xFull,ubFull,gFull,mFull,timedecomp = extForm_cheat_new(pData,disData,Ω,sN,MM,ϵ,36000,noThreads);
                lbFull = getobjectivebound(mFull);
                gapdecomp = (ubFull - lbFull)/ubFull;
            else
                global sN = sNList[Ωl];
                global MM = MMList[Ωl];
                tFull,xFull,ubFull,lbFull,timeIter,treeList,timedecomp = partSolve_BB_para_share(pData,disData,Ω,sN,MM,noThreads,5,6,5,ϵ,5,10800);
                gapdecomp = (ubFull - lbFull)/ubFull;
            end

            ubTemp = ubCalP(pData,disData,Ω,xCanBest,tCanBest,999999);
            push!(dDict[fileInd][Ωsize[Ωl]],[tFull,xFull,lbFull,ubFull,ubTemp]);
            save("test_gap.jld","lbDict",dDict,"ubDict",ubDict);
            n += 1;
            # catch
            # println("Error in Data!");
            # push!(ErrorData,(fileInd,disData));
            # end
        end
    end
end

###############
mean([dDict[fileInd][Ωsize[Ωl]][i][5] - dDict[fileInd][Ωsize[Ωl]][i][3] for i in 1:20])/(mean(ubDict[fileInd][Ωsize[Ωl]][3]))
