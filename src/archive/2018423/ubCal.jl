function ubCal(pData,disData,Ω,xhat,that)
    # calculate the upper bound of the problem given the master solution
    ubCost = that[0]*pData.p0;
    for ω in Ω
        cω = subInt(pData,disData[ω],xhat,that);
        ubCost += disData[ω].prDis*cω;
    end
    return ubCost;
end

function ubCalP(pData,disData,Ω,xhat,that)
    # parallel version of calculating the upper bound
    ubCost = that[0]*pData.p0;
    cωList = pmap(ω -> subInt(pData,disData[ω],xhat,that), Ω);
    ubCost += sum(cωList[i]*disData[Ω[i]].prDis for i in 1:length(ω));
    return ubCost;
end
