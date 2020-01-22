# extensive formulation of the multi-disruption PERT problem
global mExt = Model(with_optimizer(Gurobi.Optimizer));
global Ω = 1:3;
global HList = [10,40,70];
disData,Ωlocal = autoUGen(nameH, Hparams, nameD, dparams, 3, 1);
global dList = [disData[ω].d for ω in Ω];
global prList = [1/3, 1/3, 1/3];

function extForm(currenth, inheritData, pData, H, d, prDis, td, ωd, M = 100000, TD = 5)
    # extensive formulation could not have variable/constraint names
    # inheritData: t, x
    that = inheritData[1];
    xhat = inheritData[2];

    println("========= Disruption time $(currenth), scenario $(ωd) modeling =========");
    if td == 0
        tDict = Dict();
        xDict = Dict();
        for i in pData.II
            tDict[i] = @variable(mExt, lower_bound = 0, base_name = "t_$(td)_$(ωd)_$(i)");
            for j in pData.Ji[i]
                xDict[i,j] = @variable(mExt, lower_bound = 0, upper_bound = 1, base_name = "x_$(td)_$(ωd)_$(i)_$(j)");
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
            tDict[i] = @variable(mExt, lower_bound = 0, base_name = "t_$(td)_$(ωd)_$(i)");
            for j in pData.Ji[i]
                xDict[i,j] = @variable(mExt, lower_bound = 0, upper_bound = 1, base_name = "x_$(td)_$(ωd)_$(i)_$(j)");
                sDict[i,j] = @variable(mExt, lower_bound = 0, upper_bound = 1, base_name = "s_$(td)_$(ωd)_$(i)_$(j)");
            end
            GDict[i] = @variable(mExt, base_name = "G_$(td)_$(ωd)_$(i)", binary=true);
        end
        h = currenth + H;

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
    if td <= TD
        td += 1;
        for ω in Ω
            # sample an H
            Hnew = HList[ω];
            # sample a d
            dnew = dList[ω];
            # create new probability
            probnew = prDis*prList[ω];
            global mExt = extForm(h, [tDict,xDict], pData, Hnew, dnew, probnew, td+1, ω, M, TD);
        end
    end
    return mExt;
end
