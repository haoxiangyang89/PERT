# process of mixed first stage
@everywhere using JuMP,Gurobi,CPLEX;
@everywhere using Distributions,HDF5,JLD;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("master.jl");
@everywhere include("sub.jl");
@everywhere include("cuts.jl");
@everywhere include("detForm.jl");
@everywhere include("extForm.jl");
@everywhere include("ubCalFunc.jl");

pInputAdd = "test_14_P.csv";
kInputAdd = "test_14_K.csv";
ΩInputAdd = "test_14_Omega_full.csv";
ϕInputAdd = "test_14_Phi_full.csv";

pData = readInP(pInputAdd,kInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
data = load("testData_500.jld");
disData = data["disData"];
Ω = 1:500;
text = data["text"];
xext = data["xext"];
Gext = data["gext"];
fext = data["fext"];
GextDict = Dict();
for ω in Ω
    GextDict[ω] = Dict();
    for i in pData.II
        GextDict[ω][i] = Gext[i,ω];
    end
end

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
    γdict = Dict();
    vk = Dict();
    θInt = Dict();
    ubTemp = pData.p0*that[0];
    for ω in Ω
        Ghatω[ω] = Dict();
        for i in pData.II
            Ghatω[ω][i] = Ghat[i,ω];
        end
        πdict[ω],λdict[ω],γdict[ω],vk[ω] = subPull(pData,disData[ω],xhat,that,Ghatω[ω],400);
        ubTemp += disData[ω].prDis*subIntC(pData,disData[ω],xhat,that,400);
    end
    if ubCost > ubTemp
        ubCost = ubTemp;
        tbest = copy(that);
        xbest = copy(xhat);
    end
    ωTightCounter = 0;
    for ω in Ω
        if vk[ω] - θhat[ω] > 1e-6
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

#################################################################################
subMIntDict = Dict();
for ω in Ω
    subMIntDict[ω] = subIntGMixed(pData,disData[ω],xlb,tlb,ω,Glb);
    if θlb[ω] > subMIntDict[ω] + 1e-6
        println(θlb[ω]," ",subMIntDict[ω]);
    end
end
#################################################################################
# Print the θ^ω vs. subIntG^ω
subIntDict = Dict();
for ω in Ω
   subIntDict[ω] = subIntG(pData,disData[ω],xext,text,GextDict[ω]);
end
subMixedDict = Dict();
mMixed = solveMasterMixed(pData,disData,ωInfo,cutSet,text,xext,Gext);
solve(mMixed);
for ω in Ω
    subMixedDict[ω] = getvalue(mMixed[:θ][ω]);
end
for ω in Ω
    if subMixedDict[ω] > subIntDict[ω] + 1e-6
        println(ω," ",subMixedDict[ω]," ",subIntDict[ω]," ");
    end
end
#################################################################################

# for each i, pick an ω to separate
# Rule 1: pick the closest Hω to the current t[i]
for i in pData.II
    disMin = Inf;
    bestω = -1;
    for ω in Ω
        if !((i,ω) in ωInfo)
            if abs(disData[ω].H - tlb[i]) < disMin
                bestω = ω;
                disMin = abs(disData[ω].H - tlb[i]);
            end
        end
    end
    if bestω != -1
        push!(ωInfo,(i,bestω));
    end
end

# Rule 2: pick the mid point of fractional solution
for i in pData.II
    fracList = [ω for ω in Ω if (Glb[i,ω] < 1 - 1e-6)&(Glb[i,ω] > 1e-6)];
    if fracList != []
        startF = minimum(fracList);
        endF = maximum(fracList);
    else
        startF = -1;
        endF = -1;
    end
    println(i," ", startF," ",endF);
end
