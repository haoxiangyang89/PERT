function recoverpInfo(II,D,eff,b,K,p0,B)
    IIList = [i for i in II];
    Ddict = Dict();
    for i in 1:length(IIList)
        Ddict[IIList[i]] = D[i];
    end
    KList = [(K[k,1],K[k,2]) for k in 1:size(K)[1]];
    Succ = Dict();
    Pre = Dict();
    for i in IIList
        Succ[i] = [];
        Pre[i] = [];
    end
    for k in 1:size(K)[1]
        push!(Succ[K[k,1]],K[k,2]);
        push!(Pre[K[k,2]],K[k,1]);
    end
    bdict = Dict();
    effdict = Dict();
    Ji = Dict();
    Jmax = size(eff)[2];
    for i in 1:length(IIList)
        Ji[IIList[i]] = [j for j in 1:Jmax if eff[i,j] != 0];
        bdict[IIList[i]] = [b[i,j] for j in 1:Jmax if b[i,j] != 0];
        effdict[IIList[i]] = [eff[i,j] for j in 1:Jmax if eff[i,j] != 0];
    end
    pd = pInfo(IIList, Ji, Ddict, bdict, effdict, B, p0, KList, Pre, Succ)
    return pd;
end

function recoverdInfo(II,HOriShare,disdShare,disPrShare,Tmax)
    Ω = [];
    disData = Dict();
    for ω in 1:length(HOriShare)
        push!(Ω,ω);
        d = Dict();
        for i in 1:length(II)
            d[II[i]] = disdShare[i,ω];
        end
        disData[ω] = disInfo(HOriShare[ω],d,disPrShare[ω]);
    end
    H = Dict();
    H[0] = 0;
    # remove the duplicate ones
    counter = 1;
    for ω in Ω
        if !(disData[ω].H in values(H))
            H[counter] = disData[ω].H;
            counter += 1;
        end
    end
    H[counter] = Tmax;
    HRev = Dict();
    for hIter in keys(H)
        HRev[H[hIter]] = hIter;
    end
    return disData,H,Ω,HRev;
end

function recoverdis(IIShare,distanceShare)
    allSucc = Dict();
    distanceDict = Dict();
    for i in IIShare
        allSucc[i] = [];
    end
    for i in 1:size(distanceShare)[1]
        push!(allSucc[Int(distanceShare[i,1])],Int(distanceShare[i,2]));
        distanceDict[Int(distanceShare[i,1]),Int(distanceShare[i,2])] = distanceShare[i,3];
    end
    return allSucc,distanceDict;
end

function recoverlDict(II,lDictShare)
    lDict = Dict();
    for i in 1:length(II)
        lDict[II[i]] = lDictShare[i];
    end
    return lDict;
end

function convertDiv(divSet,divDet)
    IPPair = [(i,par) for i in pData.II for par in 1:length(divSet[i])];
    divInfoShare = SharedArray{Int,2}((length(IPPair),5));
    for ip in 1:length(IPPair)
        divInfoShare[ip,1] = IPPair[ip][1];
        divInfoShare[ip,2] = IPPair[ip][2];
        divInfoShare[ip,3] = divSet[divInfoShare[ip,1]][divInfoShare[ip,2]].startH;
        divInfoShare[ip,4] = divSet[divInfoShare[ip,1]][divInfoShare[ip,2]].endH;
        divInfoShare[ip,5] = divDet[divInfoShare[ip,1]][divInfoShare[ip,2]];
    end
    return divInfoShare;
end

function recoverDiv(divInfoShare)
    divSet = Dict();
    divDet = Dict();
    for di in 1:size(divInfoShare)[1]
        if !(divInfoShare[di,1] in keys(divSet))
            divSet[divInfoShare[di,1]] = [partType(divInfoShare[di,3],divInfoShare[di,4])];
            divDet[divInfoShare[di,1]] = [divInfoShare[di,5]];
        else
            push!(divSet[divInfoShare[di,1]],partType(divInfoShare[di,3],divInfoShare[di,4]));
            push!(divDet[divInfoShare[di,1]],divInfoShare[di,5]);
        end
    end
    return divSet,divDet;
end

function recoverCoreList(II,IJPair,textShare,xextShare,ubextShare)
    textList = [];
    xextList = [];
    ubextList = [];
    for it in 1:size(textShare)[2]
        tDict = Dict();
        for i in 1:length(II)
            tDict[II[i]] = textShare[i,it];
        end
        push!(textList,tDict);
    end
    for ijt in 1:size(xextShare)[2]
        xDict = Dict();
        for ij in 1:length(IJPair)
            xDict[IJPair[ij]] = xextShare[ij,ijt];
        end
        push!(xextList,xDict);
    end
    for iu in 1:length(ubextShare)
        push!(ubextList,ubextShare[iu]);
    end
    return textList,xextList,ubextList;
end

function runRecover(IIShare,DShare,effShare,bShare,KShare,p0BShare,HOriShare,disdShare,disPrShare,distanceShare,lDictShare)
    # recover the information and make them everywhere
    global pData = recoverpInfo(IIShare,DShare,effShare,bShare,KShare,p0BShare[1],p0BShare[2]);
    global (disData, H, Ω, HRev) = recoverdInfo(IIShare,HOriShare,disdShare,disPrShare,p0BShare[3]);
    global (allSucc,distanceDict) = recoverdis(IIShare,distanceShare);
    global lDict = recoverlDict(IIShare,lDictShare);
end
