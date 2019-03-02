function buildDistrn(nameDistr,paramDistr)
    if nameDistr == "Exponential"
        μ = paramDistr[1];
        distrObj = Exponential(μ);
    elseif nameDistr == "LogNormal"
        μ = paramDistr[1];
        σ = paramDistr[2];
        distrObj = LogNormal(μ,σ);
    elseif nameDistr == "Gamma"
        α = paramDistr[1];
        θ = paramDistr[2];
        distrObj = Gamma(α,θ);
    elseif nameDistr == "Uniform"
        la = paramDistr[1];
        ub = paramDistr[2];
        distrObj = Uniform(la,ub+1e-10);
    elseif nameDistr == "Normal"
        μ = paramDistr[1];
        σ = paramDistr[2];
        distrObj = Normal(μ,σ);
    elseif nameDistr == "PiecewiseU"
        endPts = paramDistr[1];
        mass = paramDistr[2];
        distrObj = piecewiseUniformSampler(endPts,mass);
    elseif nameDistr == "Singleton"
        endPts = [paramDistr[1],paramDistr[1]];
        mass = [1];
        distrObj = piecewiseUniformSampler(endPts,mass);
    elseif nameDistr == "Categorical"
        ps = paramDistr[1];
        xs = paramDistr[2];
        distrObj = CategoricalSamplerNew(ps,xs);
    end

    return distrObj;
end

# This is the input function that reads in the project data
function readInP(pInputAdd,kInputAdd)
    pRaw = readdlm(pInputAdd,',',header = false);
    kRaw = readdlm(kInputAdd,',',header = false);

    np,mp = size(pRaw);
    nk,tempk = size(kRaw);
    # total budget
    B = pRaw[1,1];
    # nominal scenario probability
    p0 = pRaw[1,2];

    # 0 as the last activity
    II = [0];
    Ji = Dict();
    D = Dict();
    D[0] = 0;
    Ji[0] = [];

    b = Dict();
    b[0] = Dict();
    eff = Dict();
    eff[0] = Dict();

    # read in the activity information
    for i in 2:np
        lineID = Int64(pRaw[i,1]);
        push!(II,lineID);
        D[lineID] = pRaw[i,2];

        jStart = 3;
        Ji[lineID] = [];
        b[lineID] = Dict();
        eff[lineID] = Dict();
        jCounter = 0;
        while jStart <= mp
            jCounter += 1;
            push!(Ji[lineID],jCounter);
            b[lineID][jCounter] = Float64(pRaw[i,1+2*jCounter]);
            eff[lineID][jCounter] = Float64(pRaw[i,2+2*jCounter]);
            jStart += 2;
        end
    end

    # read in the activity precedence information
    K = [];
    Pre = Dict();
    Pre[0] = [];
    Succ = Dict();
    Succ[0] = [];
    for i in II
        if i != 0
            push!(K,(i,0));
            push!(Pre[0],i);
            Pre[i] = [];
            Succ[i] = [0];
        end
    end

    for k in 1:nk
        fromI = Int64(kRaw[k,1]);
        toI = Int64(kRaw[k,2]);
        push!(K,(fromI,toI));
        push!(Pre[toI],fromI);
        push!(Succ[fromI],toI);
    end

    pData = pInfo(II,Ji,D,b,eff,B,p0,K,Pre,Succ);
    return pData;
end

# read in the scenario information from file
function readInDis(ΩInputAdd)
    ΩRaw = readdlm(ΩInputAdd,',',header = false);

    Ωn,Ωm = size(ΩRaw);
    Ω = 1:(Ωm-1);
    disData = Dict();
    for ω in Ω
        d = Dict();
        H = Float64(ΩRaw[1,ω+1]);
        pω = Float64(ΩRaw[2,ω+1]);
        for n in 3:Ωn
            d[Int64(ΩRaw[n,1])] = Float64(ΩRaw[n,ω+1]);
        end
        disData[ω] = disInfo(H,d,pω);
    end
    return disData,Ω;
end

# read in the distribution information from file
function readInUnc(ϕInputAdd)
    ϕRaw = readdlm(ϕInputAdd,',',header = false);
    ϕn,ϕm = size(ϕRaw);
    nameD = ϕRaw[1,1];
    dparams = Dict();
    for n in 2:ϕn
        dparams[Int64(ϕRaw[n,1])] = ϕRaw[n,2:ϕm];
    end
    return nameD,dparams;
end

# automatically uncertainty data generation
function autoUGen(nameH, Hparams, nameD, dparams, Ωn, totalProb)
    disData = Dict();
    Ω = 1:Ωn;
    distrH = buildDistrn(nameH,Hparams);

    distrD = Dict();
    for i in keys(dparams)
        if typeof(nameD) == Dict{Any,Any}
            distrD[i] = buildDistrn(nameD[i],dparams[i]);
        else
            distrD[i] = buildDistrn(nameD,dparams[i]);
        end
    end

    for ω in Ω
        H = round(rand(distrH),4);
        d = Dict();
        for i in keys(dparams)
            d[i] = round(rand(distrD[i]),4);
        end
        pω = totalProb/Ωn;
        disData[ω] = disInfo(H,d,pω);
    end
    return disData,Ω;
end

# automatically uncertainty data generation: stratified
function autoUGenStrata(nameH, Hparams, nameD, dparams, Ωt, Ωd, totalProb)
    disData = Dict();
    Ω = 1:(Ωt*Ωd);
    distrH = buildDistrn(nameH,Hparams);

    distrD = Dict();
    for i in keys(dparams)
        distrD[i] = buildDistrn(nameD,dparams[i]);
    end

    ω = 0;
    Ω = 1:(Ωt*Ωd);
    for ωt in 1:Ωt
        H = round(rand(distrH),4);
        for ωd in 1:Ωd
            ω += 1;
            d = Dict();
            for i in keys(dparams)
                d[i] = round(rand(distrD[i]),4);
            end
            pω = totalProb/(Ωt*Ωd);
            disData[ω] = disInfo(H,d,pω);
        end
    end
    return disData,Ω;
end

# automatically uncertainty data generation: stratified
function autoUGenStrata2(nameH, Hparams, nameD, dparams, Ωt, Ωd, totalProb)
    disData = Dict();

    distrH = buildDistrn(nameH,Hparams);

    distrD = Dict();
    for i in keys(dparams)
        distrD[i] = buildDistrn(nameD,dparams[i]);
    end

    ω = 0;
    Ω = 1:(Ωt*Ωd);
    dDict = Dict();
    for ωd in 1:Ωd
        d = Dict();
        for i in keys(dparams)
            d[i] = round(rand(distrD[i]),4);
        end
        dDict[ωd] = d;
    end
    for ωt in 1:Ωt
        H = round(rand(distrH),4);
        for ωd in 1:Ωd
            ω += 1;
            pω = totalProb/(Ωt*Ωd);
            disData[ω] = disInfo(H,dDict[ωd],pω);
        end
    end
    return disData,Ω;
end

function readH(hInputAdd)
    hRaw = readdlm(hInputAdd,',',header = false);
    nameH = hRaw[1,1];
    nh,mh = size(hRaw);
    Hparams = [];
    for i in 1:mh
        push!(Hparams,hRaw[2,i]);
    end
    return nameH,Hparams;
end

function orderdisData(disData,Ω)
    dHList = [disData[ω].H for ω in Ω];
    ωOrdered = sortperm(dHList);
    disDataNew = Dict();
    for ω in Ω
        disDataNew[ω] = disData[ωOrdered[ω]];
    end
    return disDataNew;
end

function genData(filePath,Ωsize,dataSize = 1,pName = "test_P.csv",kName = "test_K.csv",ϕName = "test_Phi.csv",hName = "test_H.csv", dOnly = 0, hOnly = 0, saveOpt = 0)
    # under the file path, look for test_P, test_K and test_Phi
    pInputAdd = joinpath(filePath,pName);
    kInputAdd = joinpath(filePath,kName);
    ϕInputAdd = joinpath(filePath,ϕName);
    hInputAdd = joinpath(filePath,hName);

    pData = readInP(pInputAdd,kInputAdd);
    nameD,dparams = readInUnc(ϕInputAdd);
    nameH,Hparams = readH(hInputAdd);

    disDataSet = [];

    for ds in 1:dataSize
        if (dOnly == 0)&(hOnly == 0)
            disData,Ω = autoUGen(nameH,Hparams,nameD,dparams,Ωsize,1 - pData.p0);
            disData = orderdisData(disData,Ω);
        elseif (dOnly != 0)&(hOnly == 0)
            distrD = Dict();
            for i in pData.II
                distrD[i] = mean(buildDistrn(nameD,dparams[i]));
            end
            disData,Ω = autoUGen(nameH,Hparams,"Singleton",distrD,Ωsize,1 - pData.p0);
            disData = orderdisData(disData,Ω);
        elseif (dOnly == 0)&(hOnly != 0)
            distrH = mean(buildDistrn(nameH,Hparams));
            disData,Ω = autoUGen("Singleton",distrH,nameD,dparams,Ωsize,1 - pData.p0);
            disData = orderdisData(disData,Ω);
        end
        push!(disDataSet,disData);
    end
    if saveOpt == 0
        return pData,disDataSet,nameD,nameH,dparams,Hparams;
    else
        save("test.jld","pData",pData,"disDataSet",disDataSet);
    end
end

function loadData(filePath,fName = "test.jld")
    # under the file path, look for test.jld (or stated specifically otherwise)
    fileAdd = joinpath(filePath,fName);
    disDataSetR = load(fileAdd);
    disDataSet = disDataSetR["disDataSet"];
    return disDataSet;
end

function genDataStrata(filePath,startaNo,sampleNo,dataSize = 1,pName = "test_P.csv",kName = "test_K.csv",ϕName = "test_Phi.csv",hName = "test_H.csv", saveOpt = 0)
    # under the file path, look for test_P, test_K and test_Phi
    pInputAdd = joinpath(filePath,pName);
    kInputAdd = joinpath(filePath,kName);
    ϕInputAdd = joinpath(filePath,ϕName);
    hInputAdd = joinpath(filePath,hName);

    pData = readInP(pInputAdd,kInputAdd);
    nameD,dparams = readInUnc(ϕInputAdd);
    nameH,Hparams = readH(hInputAdd);

    disDataSet = [];

    for ds in 1:dataSize
        disData,Ω = autoUGenStrata(nameH,Hparams,nameD,dparams,startaNo,sampleNo,1 - pData.p0);
        disData = orderdisData(disData,Ω);
        push!(disDataSet,disData);
    end
    if saveOpt == 0
        return pData,disDataSet,nameD,nameH,dparams,Hparams;
    else
        save("test.jld","pData",pData,"disDataSet",disDataSet);
    end
end
