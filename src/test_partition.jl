# test the new partition restriction method for PERT
@everywhere using JuMP,Gurobi,CPLEX;
@everywhere using Distributions,HDF5,JLD;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("master.jl");
@everywhere include("sub.jl");
@everywhere include("cuts.jl");
@everywhere include("iSolve.jl");
@everywhere include("tighten.jl");
@everywhere include("branchFunc.jl");
@everywhere include("detFunc_1l.jl");
@everywhere include("extForm_1l.jl");
@everywhere include("main_pull.jl");

pInputAdd = "test_14_P.csv";
kInputAdd = "test_14_K.csv";
ΩInputAdd = "test_14_Omega_full.csv";
ϕInputAdd = "test_14_Phi_full.csv";

pData = readInP(pInputAdd,kInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
disData,Ω = autoUGen("LogNormal",[log(35),0.5],nameD,dparams,500,1 - pData.p0);
disData = orderdisData(disData,Ω);

# create the initial partition
brInfo = precludeRel(pData,disData,Ω);
partCurrent = Dict();
partDet = Dict();
for i in pData.II
    # for each activity, break it into two partition, brInfo == 0 and brInfo == 1
    partCurrent[i] = [];
    partDet[i] = [];
    part1 = [ω for ω in Ω if brInfo[findin(pData.II,i)[1],ω] == 1];
    if part1 != []
        push!(partCurrent[i],part1);
        push!(partDet[i],1);
    end
    part0 = [ω for ω in Ω if brInfo[findin(pData.II,i)[1],ω] == 0];
    if part0 != []
        push!(partCurrent[i],part0);
        push!(partDet[i],0);
    end
end

partRev = Dict();
for i in pData.II
    partRev[i] = Dict();
    partNo = length(partCurrent[i]);
    for partIter in 1:partNo
        for item in partCurrent[i][partIter]
            partRev[i][item] = partIter;
        end
    end
end
mp = createMaster_Par(pData,disData,Ω,partCurrent,partDet);
ublb = 0;
# initialize the cutSet as empty
cutSet = Dict();
for ω in Ω
    cutSet[ω] = [];
end

stopBool = true;
while stopBool
    stopBool = false;
    mpstatus = solve(mp);
    if (mpstatus == :Optimal) && (ublb < getobjectivevalue(mp) - 1e-4)
        ublb = getobjectivevalue(mp);
        # get master solution
        that = Dict();
        xhat = Dict();
        Ghat = Dict();
        for ω in Ω
            Ghat[ω] = Dict();
            for i in pData.II
                Ghat[ω][i] = getvalue(mp[:G][i,partRev[i][ω]]);
            end
        end
        for i in pData.II
            that[i] = getvalue(mp[:t][i]);
            for j in pData.Ji[i]
                xhat[i,j] = getvalue(mp[:x][i,j]);
            end
        end

        # solve the subproblem to obtain
        πdict = Dict();
        γdict = Dict();
        λdict = Dict();
        vk = Dict();
        for ω in Ω
            πdict[ω],γdict[ω],λdict[ω],vk[ω] = subPull(pData,disData[ω],xhat,that,Ghat[ω],400);
        end
        # add the cuts
        mp,cutSet = addCuts(pData,Ω,mp,πdict,γdict,λdict,vk,that,xhat,Ghat,cutSet,partRev);
        stopBool = true;
    end
end

# obtain an upperbound
ub = getobjectivevalue(mp);
tbest = Dict();
xbest = Dict();
Gbest = Dict();
for ω in Ω
    Gbest[ω] = Dict();
    for i in pData.II
        Gbest[ω][i] = getvalue(mp[:G][i,partRev[i][ω]]);
    end
end
for i in pData.II
    tbest[i] = getvalue(mp[:t][i]);
    for j in pData.Ji[i]
        xbest[i,j] = getvalue(mp[:x][i,j]);
    end
end
