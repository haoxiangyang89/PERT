# extensive formulation of the multi-disruption PERT problem
global branchNo = 0;
function extForm(currenth, inheritData, pData, H, d, prDis, td, ωd, M = 100000, TD = 3)
    # extensive formulation could not have variable/constraint names
    # inheritData: t, x
    that = inheritData[1];
    xhat = inheritData[2];

    if td == 0
        tDict = Dict();
        xDict = Dict();
        for i in pData.II
            tDict[i] = @variable(mExt, lower_bound = 0, base_name = "t_$(branchNo)_$(i)");
            for j in pData.Ji[i]
                xDict[i,j] = @variable(mExt, lower_bound = 0, upper_bound = 1, base_name = "x_$(branchNo)_$(i)_$(j)");
            end
        end
        @constraint(mExt, sum(sum(pData.b[i][j]*xDict[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
        for k in pData.K
            @constraint(mExt, tDict[k[2]] - tDict[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*xDict[k[1],j] for j in pData.Ji[k[1]])));
        end
        for i in pData.II
            @constraint(mExt, sum(xDict[i,j] for j in pData.Ji[i]) <= 1);
        end
        @objective(mExt, Min, pData.p0*tDict[0]);
        prDis = 1 - pData.p0;
        h = currenth;
    else
        tDict = Dict();
        xDict = Dict();
        GDict = Dict();
        sDict = Dict();
        for i in pData.II
            tDict[i] = @variable(mExt, lower_bound = 0, base_name = "t_$(branchNo)_$(i)");
            for j in pData.Ji[i]
                xDict[i,j] = @variable(mExt, lower_bound = 0, upper_bound = 1, base_name = "x_$(branchNo)_$(i)_$(j)");
                sDict[i,j] = @variable(mExt, lower_bound = 0, upper_bound = 1, base_name = "s_$(branchNo)_$(i)_$(j)");
            end
            GDict[i] = @variable(mExt, base_name = "G_$(branchNo)_$(i)", binary=true);
        end
        h = currenth + H;
        println("========= Branch $(branchNo), $(td)-th disruption, disruption time $(h), scenario $(ωd) modeling =========");

        for i in pData.II
            # add indicator constraints whether t starts before h or after h
            @constraint(mExt, h - (1 - GDict[i])*M <= that[i]);
            @constraint(mExt, h + GDict[i]*M >= that[i]);
            for k in pData.Succ[i]
                # add the predecessors and the successors logic constraints
                @constraint(mExt, GDict[i] <= GDict[k]);
            end

            # add the basic sub problem constraints for the undecided activities
            @constraint(mExt, tDict[i] >= h*GDict[i]);
            @constraint(mExt, tDict[i] + GDict[i]*M >= that[i]);
            @constraint(mExt, tDict[i] - GDict[i]*M <= that[i]);

            for j in pData.Ji[i]
                @constraint(mExt, xDict[i,j] + GDict[i] >= xhat[i,j]);
                @constraint(mExt, xDict[i,j] - GDict[i] <= xhat[i,j]);

                # linearize the bilinear term of x[i,j]*G[i]
                @constraint(mExt, sDict[i,j] <= GDict[i]);
                @constraint(mExt, sDict[i,j] <= xDict[i,j]);
                @constraint(mExt, sDict[i,j] >= xDict[i,j] - 1 + GDict[i]);
            end
            @constraint(mExt, sum(xDict[i,j] for j in pData.Ji[i]) <= 1);
        end

        @constraint(mExt, sum(sum(pData.b[i][j]*xDict[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
        for k in pData.K
            @constraint(mExt, tDict[k[2]] - tDict[k[1]] >= pData.D[k[1]] + d[k[1]]*GDict[k[1]]
                - sum(pData.D[k[1]]*pData.eff[k[1]][j]*xDict[k[1],j] + d[k[1]]*pData.eff[k[1]][j]*sDict[k[1],j] for j in pData.Ji[k[1]]));
        end

        # update the objective
        objExpr = objective_function(mExt);
        objExpr += prDis*tDict[0];
        @objective(mExt, Min, objExpr);
    end

    # recursion part
    if td < TD
        for ω in Ω
            # sample an H
            Hnew = HList[ω];
            # sample a d
            dnew = dList[ω];
            # create new probability
            probnew = prDis*prList[ω];
            global branchNo += 1;
            global mExt = extForm(h, [tDict,xDict], pData, Hnew, dnew, probnew, td+1, ω, M, TD);
        end
    end
    return mExt;
end

function extForm_lin(currenth, inheritData, pData, H, d, prDis, td, ωd, M = 100000, TD = 3)
    # extensive formulation could not have variable/constraint names
    # inheritData: t, x
    that = inheritData[1];
    xhat = inheritData[2];

    if td == 0
        tDict = Dict();
        xDict = Dict();
        for i in pData.II
            tDict[i] = @variable(mExt1, lower_bound = 0, base_name = "t_$(branchNo)_$(i)");
            for j in pData.Ji[i]
                xDict[i,j] = @variable(mExt1, lower_bound = 0, upper_bound = 1, base_name = "x_$(branchNo)_$(i)_$(j)");
            end
        end
        @constraint(mExt1, sum(sum(pData.b[i][j]*xDict[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
        for k in pData.K
            @constraint(mExt1, tDict[k[2]] - tDict[k[1]] >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*xDict[k[1],j] for j in pData.Ji[k[1]])));
        end
        for i in pData.II
            @constraint(mExt1, sum(xDict[i,j] for j in pData.Ji[i]) <= 1);
        end
        @objective(mExt1, Min, pData.p0*tDict[0]);
        prDis = 1 - pData.p0;
        h = currenth;
    else
        tDict = Dict();
        xDict = Dict();
        GDict = Dict();
        sDict = Dict();
        for i in pData.II
            tDict[i] = @variable(mExt1, lower_bound = 0, base_name = "t_$(branchNo)_$(i)");
            for j in pData.Ji[i]
                xDict[i,j] = @variable(mExt1, lower_bound = 0, upper_bound = 1, base_name = "x_$(branchNo)_$(i)_$(j)");
                sDict[i,j] = @variable(mExt1, lower_bound = 0, upper_bound = 1, base_name = "s_$(branchNo)_$(i)_$(j)");
            end
            GDict[i] = @variable(mExt1, lower_bound = 0, upper_bound = 1, base_name = "G_$(branchNo)_$(i)");
        end
        h = currenth + H;
        println("========= Branch $(branchNo), $(td)-th disruption, disruption time $(h), scenario $(ωd) modeling =========");

        for i in pData.II
            # add indicator constraints whether t starts before h or after h
            @constraint(mExt1, h - (1 - GDict[i])*M <= that[i]);
            @constraint(mExt1, h + GDict[i]*M >= that[i]);
            for k in pData.Succ[i]
                # add the predecessors and the successors logic constraints
                @constraint(mExt1, GDict[i] <= GDict[k]);
            end

            # add the basic sub problem constraints for the undecided activities
            @constraint(mExt1, tDict[i] >= h*GDict[i]);
            @constraint(mExt1, tDict[i] + GDict[i]*M >= that[i]);
            @constraint(mExt1, tDict[i] - GDict[i]*M <= that[i]);

            for j in pData.Ji[i]
                @constraint(mExt1, xDict[i,j] + GDict[i] >= xhat[i,j]);
                @constraint(mExt1, xDict[i,j] - GDict[i] <= xhat[i,j]);

                # linearize the bilinear term of x[i,j]*G[i]
                @constraint(mExt1, sDict[i,j] <= GDict[i]);
                @constraint(mExt1, sDict[i,j] <= xDict[i,j]);
                @constraint(mExt1, sDict[i,j] >= xDict[i,j] - 1 + GDict[i]);
            end
            @constraint(mExt1, sum(xDict[i,j] for j in pData.Ji[i]) <= 1);
        end

        @constraint(mExt1, sum(sum(pData.b[i][j]*xDict[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
        for k in pData.K
            @constraint(mExt1, tDict[k[2]] - tDict[k[1]] >= pData.D[k[1]] + d[k[1]]*GDict[k[1]]
                - sum(pData.D[k[1]]*pData.eff[k[1]][j]*xDict[k[1],j] + d[k[1]]*pData.eff[k[1]][j]*sDict[k[1],j] for j in pData.Ji[k[1]]));
        end

        # update the objective
        objExpr = objective_function(mExt1);
        objExpr += prDis*tDict[0];
        @objective(mExt1, Min, objExpr);
    end

    # recursion part
    if td < TD
        for ω in Ω
            # sample an H
            Hnew = HList[ω];
            # sample a d
            dnew = dList[ω];
            # create new probability
            probnew = prDis*prList[ω];
            global branchNo += 1;
            global mExt1 = extForm_lin(h, [tDict,xDict], pData, Hnew, dnew, probnew, td+1, ω, M, TD);
        end
    end
    return mExt1;
end
