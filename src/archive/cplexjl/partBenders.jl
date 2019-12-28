function partBenders(cb)
    currentLB = MathProgBase.cbgetbestbound(cb);
    println("lazy,$(currentLB)");
    if currentLB <= minimum(ubCostList)
        # the callback function
        that = Dict();
        xhat = Dict();
        θhat = Dict();
        yhat = Dict();
        # obtain the solution at the current node
        for i in pData.II
            that[i] = getvalue(t[i]);
            for j in pData.Ji[i]
                xhat[i,j] = getvalue(x[i,j]);
            end
            for par in 1:length(divSet[i])
                yhat[i,par] = round(getvalue(y[i,par]));
            end
        end
        for ω in Ω
            θhat[ω] = getvalue(θ[ω]);
        end
        push!(tcoreList,that);
        push!(xcoreList,xhat);
        push!(ycoreList,yhat);

        # generate cuts
        πdict = Dict();
        λdict = Dict();
        γdict = Dict();
        vk = Dict();
        θInt = Dict();
        ubCost = minimum(ubCostList);
        ubTemp,θInt = ubCalP(pData,disData,Ω,xhat,that,Tmax1,1);
        if ubCost > ubTemp
            for i in pData.II
                tbest[i] = that[i];
                for j in pData.Ji[i]
                    xbest[i,j] = xhat[i,j];
                end
            end
        end
        push!(ubCostList,ubTemp);

        #dataList = pmap(ω -> sub_divT(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict), Ω);
        # obtain the cores
        tcore,xcore,ycore = avgCore(pData,divSet,tcoreList,xcoreList,ycoreList);
        dataList = pmap(ω -> sub_divTDualT2(pData,disData[ω],ω,that,xhat,yhat,divSet,H,lDict,tcore,xcore,ycore), Ω);
        for ω in Ω
            πdict[ω] = dataList[ω][1];
            λdict[ω] = dataList[ω][2];
            γdict[ω] = dataList[ω][3];
            vk[ω] = dataList[ω][4];
        end
        cutDual = [];
        for ω in Ω
            if vk[ω] - θhat[ω] > 1e-4*θhat[ω]
                push!(cutDual,[ω,vk[ω],πdict[ω],λdict[ω],γdict[ω]]);
                @lazyconstraint(cb, θ[ω] >= vk[ω] + sum(πdict[ω][i]*(t[i] - that[i]) for i in pData.II) +
                    sum(sum(λdict[ω][i,j]*(x[i,j] - xhat[i,j]) for j in pData.Ji[i]) for i in pData.II) +
                    sum(sum(γdict[ω][i,par]*(y[i,par] - yhat[i,par]) for par in 1:length(divSet[i])) for i in pData.II));
            end
        end
        push!(cutSet,[[that,xhat,yhat,divSet],cutDual]);
        GCurrent = [dataList[ω][5] for ω in Ω];
        push!(GList,GCurrent);
    else
        return JuMP.StopTheSolver;
    end
end
