using Distributions,HDF5,JLD;

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
        jStart = 3;
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

    if nameH == "Exponential"
        μ = Hparams[1];
        distrH = Exponential(μ);
    elseif nameH == "LogNormal"
        μ = Hparams[1];
        σ = Hparams[2];
        distrH = LogNormal(μ,σ);
    elseif nameH == "Gamma"
        α = Hparams[1];
        θ = Hparams[2];
        distrH = Gamma(α,θ);
    elseif nameH == "Uniform"
        la = Hparams[1];
        ub = Hparams[2];
        distrH = Uniform(la,ub+1e-10);
    elseif nameH == "Normal"
        μ = Hparams[1];
        σ = Hparams[2];
        distrH = Normal(μ,σ);
    end

    distrD = Dict();
    for i in keys(dparams)
        if nameD == "Exponential"
            μ = dparams[i][1];
            distrD[i] = Exponential(μ);
        elseif nameD == "LogNormal"
            μ = dparams[i][1];
            σ = dparams[i][2];
            distrD[i] = LogNormal(μ,σ);
        elseif nameD == "Gamma"
            α = dparams[i][1];
            θ = dparams[i][2];
            distrD[i] = Gamma(α,θ);
        elseif nameD == "Uniform"
            la = dparams[i][1];
            ub = dparams[i][2];
            distrD[i] = Uniform(la,ub+1e-10);
        elseif nameD == "Normal"
            μ = dparams[i][1];
            σ = dparams[i][2];
            distrD[i] = Normal(μ,σ);
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

    if nameH == "Exponential"
        μ = Hparams[1];
        distrH = Exponential(μ);
    elseif nameH == "LogNormal"
        μ = Hparams[1];
        σ = Hparams[2];
        distrH = LogNormal(μ,σ);
    elseif nameH == "Gamma"
        α = Hparams[1];
        θ = Hparams[2];
        distrH = Gamma(α,θ);
    elseif nameH == "Uniform"
        la = Hparams[1];
        ub = Hparams[2];
        distrH = Uniform(la,ub+1e-10);
    elseif nameH == "Normal"
        μ = Hparams[1];
        σ = Hparams[2];
        distrH = Normal(μ,σ);
    end

    distrD = Dict();
    for i in keys(dparams)
        if nameD == "Exponential"
            μ = dparams[i][1];
            distrD[i] = Exponential(μ);
        elseif nameD == "LogNormal"
            μ = dparams[i][1];
            σ = dparams[i][2];
            distrD[i] = LogNormal(μ,σ);
        elseif nameD == "Gamma"
            α = dparams[i][1];
            θ = dparams[i][2];
            distrD[i] = Gamma(α,θ);
        elseif nameD == "Uniform"
            la = dparams[i][1];
            ub = dparams[i][2];
            distrD[i] = Uniform(la,ub+1e-10);
        elseif nameD == "Normal"
            μ = dparams[i][1];
            σ = dparams[i][2];
            distrD[i] = Normal(μ,σ);
        end
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

    if nameH == "Exponential"
        μ = Hparams[1];
        distrH = Exponential(μ);
    elseif nameH == "LogNormal"
        μ = Hparams[1];
        σ = Hparams[2];
        distrH = LogNormal(μ,σ);
    elseif nameH == "Gamma"
        α = Hparams[1];
        θ = Hparams[2];
        distrH = Gamma(α,θ);
    elseif nameH == "Uniform"
        la = Hparams[1];
        ub = Hparams[2];
        distrH = Uniform(la,ub+1e-10);
    elseif nameH == "Normal"
        μ = Hparams[1];
        σ = Hparams[2];
        distrH = Normal(μ,σ);
    end

    distrD = Dict();
    for i in keys(dparams)
        if nameD == "Exponential"
            μ = dparams[i][1];
            distrD[i] = Exponential(μ);
        elseif nameD == "LogNormal"
            μ = dparams[i][1];
            σ = dparams[i][2];
            distrD[i] = LogNormal(μ,σ);
        elseif nameD == "Gamma"
            α = dparams[i][1];
            θ = dparams[i][2];
            distrD[i] = Gamma(α,θ);
        elseif nameD == "Uniform"
            la = dparams[i][1];
            ub = dparams[i][2];
            distrD[i] = Uniform(la,ub+1e-10);
        elseif nameD == "Normal"
            μ = dparams[i][1];
            σ = dparams[i][2];
            distrD[i] = Normal(μ,σ);
        end
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
