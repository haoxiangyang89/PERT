# Two bus case with piecewise constant f(H)
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("master.jl");
@everywhere include("sub.jl");
@everywhere include("cuts.jl");
@everywhere include("detForm.jl");
@everywhere include("extForm.jl");
@everywhere include("ubCalFunc.jl");

pInputAdd = "test_2_P_PieceU.csv";
kInputAdd = "test_2_K.csv";
ϕInputAdd = "test_2_Phi_PieceU.csv";

pData = readInP(pInputAdd,kInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
disData,Ω = autoUGen("PiecewiseU",[[0.1,1.5,1.7],[0.2,0.8]],nameD,dparams,5000,1 - pData.p0);
disData = orderdisData(disData,Ω);
@time text,xext,fext,gext,mp = extForm_cheat(pData,disData,Ω);

################################################################################
tdet,xdet,fdet = detBuild(pData);
ubdet = ubCal(pData,disData,Ω,xdet,tdet);
brInfo = precludeRel(pData,disData,Ω,ubdet);
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
    partn1 = [ω for ω in Ω if brInfo[findin(pData.II,i)[1],ω] == -1];
    if partn1 != []
        push!(partCurrent[i],partn1);
        push!(partDet[i],-1);
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

ωInfo = [];
cutSet = Dict();
for ω in Ω
    cutSet[ω] = [];
end
ubCost = Inf;
lbCost = -Inf;

keepIter = true;
tlb = Dict();
xlb = Dict();
Glb = Dict();
θlb = Dict();
xbest = Dict();
tbest = Dict();
ωWrong = [];
while keepIter
    mp = createMaster_MixedL(pData,disData,Ω,ωInfo,cutSet,partCurrent,partDet,400);
    solve(mp);
    # obtain the solution
    that = Dict();
    xhat = Dict();
    Ghat = Dict();
    θhat = Dict();
    Ghatω = Dict();
    for i in pData.II
        that[i] = getvalue(mp[:t][i]);
        for j in pData.Ji[i]
            xhat[i,j] = getvalue(mp[:x][i,j]);
        end
        for ω in Ω
            Ghat[i,ω] = getvalue(mp[:G][i,ω]);
        end
    end
    for ω in Ω
        θhat[ω] = getvalue(mp[:θ][ω]);
    end
    lbCost = getobjectivevalue(mp);
    # generate cuts
    lbPrev = lbCost;
    πdict = Dict();
    λdict = Dict();
    πdict1 = Dict();
    λdict1 = Dict();
    γdict = Dict();
    vk = Dict();
    vk1 = Dict();
    θInt = Dict();
    for ω in Ω
        Ghatω[ω] = Dict();
        for i in pData.II
            Ghatω[ω][i] = Ghat[i,ω];
        end
    end
    # for ω in Ω
    #     πdict[ω],λdict[ω],γdict[ω],vk[ω] = subPull(pData,disData[ω],xhat,that,Ghatω[ω],400);
    # end
    ubTemp = ubCal(pData,disData,Ω,xhat,that);
    dataList = pmap(ω -> subPull(pData,disData[ω],xhat,that,Ghatω[ω],400), Ω);
    for ω in Ω
        πdict[ω] = dataList[ω][1];
        λdict[ω] = dataList[ω][2];
        γdict[ω] = dataList[ω][3];
        vk[ω] = dataList[ω][4];
    end
        #πdict1[ω],λdict1[ω],vk1[ω] = subLag(pData,disData[ω],xhat,that,Ghatω[ω],γdict[ω],400);
    if ubCost > ubTemp
        ubCost = ubTemp;
        tbest = copy(that);
        xbest = copy(xhat);
    end
    ωTightCounter = 0;
    for ω in Ω
        if vk[ω] - θhat[ω] > 1e-5
            push!(cutSet[ω],(πdict[ω],λdict[ω],γdict[ω],vk[ω],that,xhat,Ghatω[ω]));
        else
            ωTightCounter += 1;
        end
    end
    if ωTightCounter == length(Ω)
        keepIter = false;
        for i in pData.II
            tlb[i] = that[i];
            for j in pData.Ji[i]
                xlb[i,j] = xhat[i,j];
            end
            for ω in Ω
                Glb[i,ω] = Ghat[i,ω];
            end
        end
        for ω in Ω
            θlb[ω] = θhat[ω];
        end
    end
end
################################################################################
t2Improve = [];
for ω in Ω
    if (Glb[2,ω] < 1 - 1e-7)&(Glb[2,ω] > 1e-7)
        lbProxy = masterMaxLP(pData,disData,Ω,ωInfo,Glb,[(i,ω)],cutSet,partCurrent,partDet,400);
        push!(t2Improve,lbProxy);
    end
end
ωSet = [ω for ω in Ω if (Glb[2,ω] < 1 - 1e-7)&(Glb[2,ω] > 1e-7)];
HSet = [disData[ω].H for ω in ωSet];

push!(ωInfo,(2,ωSet[indmax(t2Improve)]));
################################################################################
