# serial example maker
function serialMaker1(n,k,D0,d0,e0,p0,B,b0,nameH,Hparams,Ωsize)
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
    for i in 1:n
        push!(II,i);
        D[i] = k^(n-i)*D0;

        Ji[i] = [1];
        b[i] = Dict();
        b[i][1] = b0;
        eff[i] = Dict();
        eff[i][1] = e0;
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

    for i in 1:(n - 1)
        fromI = i;
        toI = i + 1;
        push!(K,(fromI,toI));
        push!(Pre[toI],fromI);
        push!(Succ[fromI],toI);
    end

    pData = pInfo(II,Ji,D,b,eff,B,p0,K,Pre,Succ);

    distrD = Dict();
    for i in 1:n
        distrD[i] = d0*k^(i - 1);
    end
    distrD[0] = 0;
    disData = Dict();
    #H = mean(buildDistrn(nameH,Hparams))
    #disData[1] = disInfo(H,distrD,1 - pData.p0)
    disData,Ω = autoUGen(nameH,Hparams,"Singleton",distrD,Ωsize,1 - pData.p0);
    disData = orderdisData(disData,Ω);

    return pData,disData,Ω;
end

function serialMaker2(n,k,D0,d0,e0,p0,B,b0,nameH,Hparams,Ωsize)
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
    for i in 1:n
        push!(II,i);
        D[i] = D0*k^(i - 1);

        Ji[i] = [1];
        b[i] = Dict();
        b[i][1] = b0;
        eff[i] = Dict();
        eff[i][1] = e0;
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

    for i in 1:(n - 1)
        fromI = i;
        toI = i + 1;
        push!(K,(fromI,toI));
        push!(Pre[toI],fromI);
        push!(Succ[fromI],toI);
    end

    pData = pInfo(II,Ji,D,b,eff,B,p0,K,Pre,Succ);

    distrD = Dict();
    for i in 1:n
        distrD[i] = d0*k^(i - 1);
    end
    distrD[0] = 0;
    disData = Dict();
    #H = mean(buildDistrn(nameH,Hparams))
    #disData[1] = disInfo(H,distrD,1 - pData.p0)
    disData,Ω = autoUGen(nameH,Hparams,"Singleton",distrD,Ωsize,1 - pData.p0);
    disData = orderdisData(disData,Ω);

    return pData,disData,Ω;
end

function specialMaker(k,k1,D0,d0,e0,p0,pd,B,b0,nameH,Hparams,Ωsize,distrName,distrD)
    # make the special case where Full << exp/HOnly/dOnly/det
    II = [0,1,2,3,4];
    for i in 1:2
        push!(II,4+i);
    end

    Pre = Dict();
    Succ = Dict();
    for i in II
        Pre[i] = [];
        Succ[i] = [];
    end

    K = [(1,2),(3,4)];
    push!(Pre[2],1);
    push!(Succ[1],2);
    push!(Pre[4],3);
    push!(Succ[3],4);
    for i in 1:2
        push!(K,(4,4+i));
        push!(Pre[4+i],4);
        push!(Succ[4],4+i);
    end

    for i in II
        if i != 0
            push!(K,(i,0));
            push!(Pre[0],i);
            push!(Succ[i],0);
        end
    end

    Ji = Dict();
    D = Dict();
    D[0] = 0;
    Ji[0] = [];

    b = Dict();
    b[0] = Dict();
    eff = Dict();
    eff[0] = Dict();

    # set up the activity duration, crashing budget and effect
    D[1] = k*D0;
    D[2] = D0;
    D[3] = k1*D0;
    D[4] = D0;
    for i in 1:2
        D[4+i] = D0;
    end
    for i in II
        Ji[i] = [1];
        b[i] = Dict();
        b[i][1] = b0;
        eff[i] = Dict();
        eff[i][1] = e0;
    end
    pData = pInfo(II,Ji,D,b,eff,B,p0,K,Pre,Succ);

    # set up the activity disruption information
    disData = Dict();
    disData,Ω = autoUGen(nameH,Hparams,distrName,distrD,Ωsize,1 - pData.p0);
    disData = orderdisData(disData,Ω);

    return pData,disData,Ω;
end

# experiments file generator
function experimentgen(filePath,dataSize,Ωsize,Hskewness,μp,Dmag,Doption = "r",pName = "test_P.csv",kName = "test_K.csv")
    # read in the deterministic network
    pInputAdd = joinpath(filePath,pName);
    kInputAdd = joinpath(filePath,kName);

    pData = readInP(pInputAdd,kInputAdd);
    nameD = "Exponential";
    nameH = "PiecewiseU";
    tdet,xdet,fdet = detBuild(pData);
    TDet = tdet[0];

    # determine the H and d distributions
    #Hparams = [[0,(1 - Hskewness)*TDet,TDet],[Hskewness,1 - Hskewness]];
    Hμ = μp*TDet;
    Hparams = [[0,Hμ,Hμ/(1 - Hskewness)],[Hskewness,1-Hskewness]];
    dparams = Dict();
    if Doption == "a"
        for i in pData.II
            if i != 0
                dparams[i] = Dmag*pData.D[i];
            end
        end
    elseif Doption == "d"
        dInd = sortperm([pData.D[i] for i in pData.II if i != 0]);
        for i in 1:length(dInd)
            dparams[dInd[i]] = Dmag*pData.D[length(dInd) - i + 1];
        end
    else
        randInd = randperm(length(pData.II) - 1);       # excluding 0
        for i in pData.II
            if i != 0
                dparams[i] = Dmag*pData.D[randInd[i]];
            end
        end
    end

    disDataSet = [];
    # generate the data and save to filePath
    for ds in 1:dataSize
        disData,Ω = autoUGen(nameH,Hparams,nameD,dparams,Ωsize,1 - pData.p0);
        disData = orderdisData(disData,Ω);
        push!(disDataSet,disData);
    end

    return pData,disDataSet;
end
