function partitionSolve(pData,disData,ϵ = 0.01,tightenBool = 0)
    # process to solve the PERT problem
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

    # start with an upper bound based on the deterministic solution
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

    while (ubCost - lbCost)/ubCost > ϵ
        keepIter = true;
        tlb = Dict();
        xlb = Dict();
        θlb = Dict();
        ylb = Dict();
        dataList = [];
        mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax);
        # if perform the bound tightening process
        if tightenBool == 1
            mpTemp = copy(mp);
            ubInfo,lbInfo = obtainBds(pData,disData,Ω,mpTemp,ubCost);
            mp = updateMaster(mp,ubInfo,lbInfo);
            divSet,divDet = revisePar(pData,disData,divSet,divDet,ubInfo,lbInfo);
            mp = createMaster_Div(pData,disData,Ω,divSet,divDet,cutSet,Tmax);
        end
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
            dataList = pmap(ω -> sub_divT(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict), Ω);
            for ω in Ω
                πdict[ω] = dataList[ω][1];
                λdict[ω] = dataList[ω][2];
                γdict[ω] = dataList[ω][3];
                vk[ω] = dataList[ω][4];
            end
            ωTightCounter = 0;
            cutDual = Dict();
            for ω in Ω
                if vk[ω] - θhat[ω] > 1e-5
                    cutDual[ω] = [vk[ω],πdict[ω],λdict[ω],γdict[ω]];
                    mp = addtxyCut(pData,ω,mp,πdict[ω],λdict[ω],γdict[ω],vk[ω],that,xhat,yhat,divSet);
                else
                    cutDual[ω] = [];
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
                    for par in 1:length(divSet[i])
                        ylb[i,par] = yhat[i,par];
                    end
                end
                for ω in Ω
                    θlb[ω] = θhat[ω];
                end
            else
                push!(cutSet,[[that,xhat,yhat,divSet],cutDual]);
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
    end

    return tbest,xbest,lbCost,ubCost;
end
