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
lbPrev = -Inf;
ωWrong = [];
while keepIter
    mp = createMaster_Mixed(pData,disData,Ω,ωInfo,cutSet,400);
    solve(mp);
    # obtain the solution
    that = Dict();
    xhat = Dict();
    Ghat = Dict();
    θhat = Dict();
    for i in pData.II
        that[i] = getvalue(mp[:t][i]);
        for j in pData.Ji[i]
            xhat[i,j] = getvalue(mp[:x][i,j]);
        end
        for ω in Ω
            if (i,ω) in ωInfo
                Ghat[i,ω] = getvalue(mp[:G][i,ω]);
            end
        end
    end
    for ω in Ω
        θhat[ω] = getvalue(mp[:θ][ω]);
    end
    lbCost = getobjectivevalue(mp);
    if lbPrev < lbCost - 1e-5
        lbPrev = lbCost;
        πdict = Dict();
        λdict = Dict();
        γdict = Dict();
        vk = Dict();
        θInt = Dict();
        ubTemp = pData.p0*that[0];
        for ω in Ω
            πdict[ω],λdict[ω],γdict[ω],vk[ω] = subMixed(pData,disData[ω],xhat,that,Ghat,ω,Ω,400);
            push!(cutSet[ω],(πdict[ω],λdict[ω],γdict[ω],vk[ω],that,xhat,Ghat));
            ubTemp += disData[ω].prDis*subIntGMixed(pData,disData[ω],xhat,that,ω,Ghat);
        end
        if ubCost > ubTemp
            ubCost = ubTemp;
            tbest = copy(that);
            xbest = copy(xhat);
        end
    else
        keepIter = false;
        for i in pData.II
            tlb[i] = that[i];
            for j in pData.Ji[i]
                xlb[i,j] = xhat[i,j];
            end
            for ω in Ω
                if (i,ω) in ωInfo
                    Glb[i,ω] = Ghat[i,ω];
                end
            end
        end
        for ω in Ω
            θlb[ω] = θhat[ω];
        end
    end
end
lbPrev = lbCost;

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
