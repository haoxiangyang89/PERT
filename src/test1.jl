addprocs(30);
global noThreads = 30;
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;
# test sbb
@everywhere include("header.jl");

pathList = ["/home/haoxiang/PERT_tests/current/11/",
            "/home/haoxiang/PERT_tests/current/14/",
            "/home/haoxiang/PERT_tests/current/19/"];
#            ,"/home/haoxiang/PERT_tests/current/35/",
#            "/home/haoxiang/PERT_tests/current/55/",
#            "/home/haoxiang/PERT_tests/current/75/"];
HsList = 0.1:0.1:0.9;
for fileInd in 1:length(pathList)
    filePath = pathList[fileInd];
    dDict = Dict();
    for Hskewness in HsList
        dataRaw = load(filePath*"test1_$(fileInd)_$(Hskewness).jld");
        pData = dataRaw["pData"];
        disDataSet = dataRaw["disDataSet"];
        dDict[Hskewness] = Dict();
        for bNo in 1:20
            disData = disDataSet[bNo];
            Ω = 1:length(disData);
            # deterministic
            tdet,xdet,fdet = detBuild(pData);
            ubdet = ubCalP(pData,disData,Ω,xdet,tdet,9999999);

            # expected
            eH = mean([disData[ω].H for ω in Ω]);
            ed = Dict();
            for i in pData.II
                if i == 0
                    ed[i] = 0;
                else
                    ed[i] = mean([disData[ω].d[i] for ω in Ω]);
                end
            end
            texp,xexp,fexp,Gexp,mexp = expModel(pData,eH,ed);
            ubexp = ubCalP(pData,disData,Ω,xexp,texp,9999999);

            # dOnly
            disData1 = deepcopy(disData);
            for ω in Ω
                disData[ω].H = eH;
            end
            tdOnly,xdOnly,fdOnly,gdOnly,mdOnly = extForm_cheat(pData,disData,Ω,1e-4,999999);
            disData = deepcopy(disData1);
            ubdOnly = ubCalP(pData,disData,Ω,xdOnly,tdOnly,999999);

            # hOnly
            for ω in Ω
                for i in pData.II
                    if i != 0
                        disData[ω].d[i] = ed[i];
                        if disData[ω].d[i] < 1e-4
                            disData[ω].d[i] = 0;
                        end
                    end
                end
            end
            tHOnly,xHOnly,fHOnly,gHOnly,mHOnly = extForm_cheat(pData,disData,Ω,1e-4,9999999);
            disData = deepcopy(disData1);
            ubHOnly = ubCalP(pData,disData,Ω,xHOnly,tHOnly,9999999);

            # full
            tFull,xFull,ubFull,gFull,mFull = extForm_cheat(pData,disData,Ω,1e-4,9999999);
            dDict[Hskewness][bNo] = [tdet,xdet,fdet,ubdet,
                                texp,xexp,fexp,ubexp,
                                tdOnly,xdOnly,fdOnly,ubdOnly,
                                tHOnly,xHOnly,fHOnly,ubHOnly,
                                tFull,xFull,ubFull];
        end
    end
    save(filePath*"test1_data.jld","dDict",dDict);
end
