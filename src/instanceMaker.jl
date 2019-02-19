# serial example maker
function serialMaker(n,k,D0,d0,e0,p0,B,nameH,Hparams,Ωsize)
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
        b[i][1] = 1;
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
    H = mean(buildDistrn(nameH,Hparams))
    #disData[1] = disInfo(H,distrD,1 - pData.p0)
    disData,Ω = autoUGen(nameH,Hparams,"Singleton",distrD,Ωsize,1 - pData.p0);

    return pData,disData;
end
