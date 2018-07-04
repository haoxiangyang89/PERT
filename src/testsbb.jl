addprocs(20);
# test sbb
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
@everywhere include("part_tight.jl");

pInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_P.csv";
kInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_K.csv";
ϕInputAdd = "/home/haoxiang/PERT_tests/14_ExponentialD_LogNormalH/test_14_Phi_full.csv";

pData = readInP(pInputAdd,kInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
disData,Ω = autoUGen("LogNormal",[log(35),0.5],nameD,dparams,500,1 - pData.p0);
disData = orderdisData(disData,Ω);

Tmax = disData[length(Ω)].H + longestPath(pData)[0];
pdData = deepcopy(pData);
for i in pData.II
    if i != 0
        pdData.D[i] = pData.D[i] + maximum([disData[ω].d[i] for ω in Ω])
    else
        pdData.D[i] = pData.D[i];
    end
end
lDict = longestPath(pdData);
for i in pData.II
    lDict[i] += disData[length(Ω)].H;
end
Tmax1 = lDict[0];
tdet,xdet,fdet = detBuild(pData);
ubdet = ubCal(pData,disData,Ω,xdet,tdet,Tmax1);
brInfo = precludeRel(pData,disData,Ω,ubdet);

H = Dict();
H[0] = 0;
H[length(Ω)+1] = Tmax;
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
ubCost = ubdet;
lbCost = -Inf;

keepIter = true;
tlb = Dict();
xlb = Dict();
θlb = Dict();
ylb = Dict();
dataList = [];
mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax);
mpTemp = copy(mp);
ubInfo,lbInfo = obtainBds(pData,disData,Ω,mpTemp,ubCost);
# divSet,divDet = revisePar(pData,disData,divSet,divDet,ubInfo,lbInfo);
# mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax);
mp = updateMaster(mp,ubInfo,lbInfo);
while keepIter
    solve(mp);
    # obtain the solution
    that = Dict();
    xhat = Dict();
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
    ubTemp,θInt = ubCalP(pData,disData,Ω,xhat,that,Tmax1,1);
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
    # dataList = pmap(ω -> sub_div(pData,disData[ω],ω,that,xhat,yhat,divSet,ubInfo,lbInfo,Tmax1), Ω);
    dataList = pmap(ω -> sub_divT(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict), Ω);
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
        if vk[ω] - θhat[ω] > 1e-5
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
GFrac = Dict();
for i in pData.II
    GFraciList = [ω for ω in Ω if (dataList[ω][5][i] < 1 - 1e-6)&(dataList[ω][5][i] > 1e-6)];
    if GFraciList != []
        GFrac[i] = [minimum(GFraciList),maximum(GFraciList)];
    else
        GFrac[i] = [];
    end
end
# create new partition
newPartition = [];
for i in pData.II
    if GFrac[i] != []
        newItem = (i,Int(floor((GFrac[i][1] + GFrac[i][2])/2)));
        push!(newPartition,newItem);
    end
end
divSet,divDet = splitPar(divSet,divDet,newPartition);
