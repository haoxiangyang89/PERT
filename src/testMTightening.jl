# test for branch and bound tightening + big M tightening
# Two bus case with piecewise constant f(H)
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("master.jl");
@everywhere include("sub.jl");
@everywhere include("cuts.jl");
@everywhere include("iSolve.jl");
@everywhere include("tighten.jl");
@everywhere include("branchFunc.jl");
@everywhere include("detForm.jl");
@everywhere include("extForm.jl");
@everywhere include("ubCalFunc.jl");
@everywhere include("tighten.jl");
@everywhere include("partition_LP.jl");
@everywhere include("partition_LR.jl");

pInputAdd = "test_2_P_PieceU.csv";
kInputAdd = "test_2_K.csv";
ϕInputAdd = "test_2_Phi_PieceU.csv";

pData = readInP(pInputAdd,kInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
disData,Ω = autoUGen("PiecewiseU",[[0.1,1.5,1.7],[0.2,0.8]],nameD,dparams,500,1 - pData.p0);
disData = orderdisData(disData,Ω);
@time text,xext,fext,gext,mp = extForm_cheat(pData,disData,Ω);

################################################################################
tdet,xdet,fdet = detBuild(pData);
ubdet = ubCal(pData,disData,Ω,xdet,tdet);
brInfo = precludeRel(pData,disData,Ω,ubdet);

################################################################################
# initialization
mp = createMaster(pData,disData,Ω);
keepIter = true;
tlb = Dict();
xlb = Dict();
θlb = Dict();
xbest = Dict();
tbest = Dict();
ubCost = ubdet;
mpTemp = copy(mp);
ubInfo,lbInfo = obtainBds(pData,disData,Ω,mpTemp,ubdet);
mp = updateMaster(mp,ubInfo,lbInfo);
# Benders decomposition
while keepIter
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
    end
    for ω in Ω
        θhat[ω] = getvalue(mp[:θ][ω]);
    end
    lbCost = getobjectivevalue(mp);
    # generate cuts
    lbPrev = lbCost;
    πdict = Dict();
    λdict = Dict();
    vk = Dict();
    θInt = Dict();
    ubTemp = ubCal(pData,disData,Ω,xhat,that);
    if ubCost > ubTemp
        ubCost = ubTemp;
        tbest = copy(that);
        xbest = copy(xhat);
    end
    # dataList = Dict();
    # for ω in Ω
    #     dataList[ω] = subTight(pData,disData[ω],xhat,that,ubInfo,lbInfo);
    #     println(ω);
    # end
    dataList = pmap(ω -> subTight(pData,disData[ω],xhat,that,ubInfo,lbInfo), Ω);
    for ω in Ω
        πdict[ω] = dataList[ω][1];
        λdict[ω] = dataList[ω][2];
        vk[ω] = dataList[ω][3];
    end
        #πdict1[ω],λdict1[ω],vk1[ω] = subLag(pData,disData[ω],xhat,that,Ghatω[ω],γdict[ω],400);
    ωTightCounter = 0;
    for ω in Ω
        if vk[ω] - θhat[ω] > 1e-5
            @constraint(mp,mp[:θ][ω] >= vk[ω] + sum(πdict[ω][i]*(mp[:t][i] - that[i]) for i in pData.II) +
                sum(sum(λdict[ω][i,j]*(mp[:x][i,j] - xhat[i,j]) for j in pData.Ji[i]) for i in pData.II));
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
        end
        for ω in Ω
            θlb[ω] = θhat[ω];
        end
    end
end

################################################################################
# spatial b&b process
H = Dict();
H[0] = 0;
H[length(Ω)+1] = 10;
for ω in Ω
    H[ω] = disData[ω].H;
end

# initialize cutSet and divSet
cutSet = [];
divSet = Dict();
divDet = Dict();
for i in pData.II
    set1 = [ω for ω in Ω if brInfo[findfirst(pData.II,i),ω] == 1];
    setn1 = [ω for ω in Ω if brInfo[findfirst(pData.II,i),ω] == -1];

    if set1 != []
        set1t = partType(0,maximum(set1));
        if setn1 != []
            setn1t = partType(minimum(setn1),length(Ω) + 1);
            set0t = partType(maximum(set1),minimum(setn1));
            divSet[i] = [set1t,set0t,setn1t];
            divDet[i] = [1,0,-1];
        else
            set0t = partType(maximum(set1),length(Ω) + 1);
            divSet[i] = [set1t,set0t];
            divDet[i] = [1,0];
        end
    else
        if setn1 != []
            setn1t = partType(minimum(setn1),length(Ω) + 1);
            set0t = partType(0,minimum(setn1));
            divSet[i] = [set0t,setn1t];
            divDet[i] = [0,-1];
        else
            set0t = partType(0,length(Ω) + 1);
            divSet[i] = [set0t];
            divDet[i] = [0];
        end
    end
end
xbest = Dict();
tbest = Dict();

keepIter = true;
tlb = Dict();
xlb = Dict();
θlb = Dict();
ylb = Dict();
mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax);
while keepIter
    solve(mp);
    # obtain the solution
    that = Dict();
    xhat = Dict();
    Ghat = Dict();
    θhat = Dict();
    yhat = Dict();
    for i in pData.II
        that[i] = getvalue(mp[:t][i]);
        for j in pData.Ji[i]
            xhat[i,j] = getvalue(mp[:x][i,j]);
        end
        for par in 1:length(divSet[i])
            yhat[i,par] = getvalue(mp[:y][i,par]);
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
    ubTemp = ubCal(pData,disData,Ω,xhat,that);
    if ubCost > ubTemp
        ubCost = ubTemp;
        tbest = copy(that);
        xbest = copy(xhat);
    end
    # dataList = Dict();
    # for ω in Ω
    #     dataList[ω] = sub_div(pData,disData[ω],ω,that,xhat,yhat,divSet,100,1);
    #     println(ω);
    # end
    dataList = pmap(ω -> sub_div(pData,disData[ω],ω,that,xhat,yhat,divSet,100), Ω);
    for ω in Ω
        πdict[ω] = dataList[ω][1];
        λdict[ω] = dataList[ω][2];
        γdict[ω] = dataList[ω][3];
        vk[ω] = dataList[ω][4];
    end
        #πdict1[ω],λdict1[ω],vk1[ω] = subLag(pData,disData[ω],xhat,that,Ghatω[ω],γdict[ω],400);
    ωTightCounter = 0;
    cutDual = Dict();
    for ω in Ω
        if vk[ω] - θhat[ω] > 1e-6
            cutDual[ω] = [vk[ω],πdict[ω],λdict[ω],γdict[ω]];
            mp = addtxCut(pData,ω,mp,πdict,λdict,γdict,vk,that,xhat,yhat,divSet);
        else
            cutDual[ω] = [];
            ωTightCounter += 1;
        end
    end
    push!(cutSet,[[that,xhat,yhat,divSet],cutDual]);
    if ωTightCounter == length(Ω)
        keepIter = false;
        for i in pData.II
            tlb[i] = that[i];
            for j in pData.Ji[i]
                xlb[i,j] = xhat[i,j];
            end
            for par in 1:length(divSet[i])
                ylb[i,par] = yhat[i,par];
            end
        end
        for ω in Ω
            θlb[ω] = θhat[ω];
        end
    end
end

# need to come up with a rule to partition: gradient descent like binary search
# check θInt vs. θhat: why the lower bound and the upper bound do not converge quickly --->
# use the sub problem solution G to learn the b&b
# also need to think up a way to tightening the cuts for each partition
divSet,divDet = splitPar(divSet,divDet,[(2,102)]);
