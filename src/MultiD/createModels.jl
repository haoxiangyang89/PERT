# add state variables
function add_state_variables(sp, pData)
    @variable(sp, 0 <= t[i in pData.II] <= 100000, SDDP.State, initial_value = 0);
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1, SDDP.State, initial_value = 0);
    @variable(sp, 0 <= h <= 50000, SDDP.State, initial_value = 0);
end


# create stage models
function create_first_stage(sp, pData)
    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*sp[:x][i,j].out for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, durationConstr[k in pData.K], sp[:t][k[2]].out - sp[:t][k[1]].out >= pData.D[k[1]]*(1-sum(pData.eff[k[1]][j]*sp[:x][k[1],j].out for j in pData.Ji[k[1]])));
    @constraint(sp, xConstr[i in pData.II], sum(sp[:x][i,j].out for j in pData.Ji[i]) <= 1);
    @constraint(sp, inih, sp[:h].out == 0);
    @stageobjective(sp, pData.p0*sp[:t][0].out);
end

function create_n_stage(sp, n, pData, dDn, M)
    @variable(sp, G[i in pData.II], Bin);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);
    # add disruption timing and magnitude as variables
    @variables(sp, begin
        H
        d[i in pData.II]
        gh[i in pData.II]
    end);
    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II],sp[:h].out - (1 - G[i])*M <= sp[:t][i].in);
    @constraint(sp, GCons[i in pData.II],sp[:h].out + G[i]*M >= sp[:t][i].in);

    # add the predecessors and the successors logic constraints
    @constraint(sp, GSuccessors[i in pData.II, k in pData.Succ[i]], G[i] <= G[k]);

    # add the linearization constraints for G*h
    @constraint(sp, ghCons1[i in pData.II], gh[i] <= G[i]*pData.hMax);
    @constraint(sp, ghCons2[i in pData.II], gh[i] <= sp[:h].out);
    @constraint(sp, ghCons3[i in pData.II], gh[i] >= G[i] + sp[:h].out - 1);

    # add h update dynamics
    @constraint(sp, hCons, sp[:h].out == sp[:h].in + H);

    # add the basic sub problem constraints for the undecided activities
    @constraint(sp, tGbound1[i in pData.II],sp[:t][i].out >= gh[i]);
    @constraint(sp, tFnAnt1[i in pData.II],sp[:t][i].out + G[i]*M >= sp[:t][i].in);
    @constraint(sp, tFnAnt2[i in pData.II],sp[:t][i].out - G[i]*M <= sp[:t][i].in);
    @constraint(sp, xFnAnt1[i in pData.II, j in pData.Ji[i]],sp[:x][i,j].out + G[i] >= sp[:x][i,j].in);
    @constraint(sp, xFnAnt2[i in pData.II, j in pData.Ji[i]],sp[:x][i,j].out - G[i] <= sp[:x][i,j].in);

    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= sp[:x][i,j].out);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= sp[:x][i,j].out - 1 + G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*sp[:x][i,j].out for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(sp[:x][i,j].out for j in pData.Ji[i]) <= 1);
    @constraint(sp, durationConstr[k in pData.K], sp[:t][k[2]].out - sp[:t][k[1]].out >= pData.D[k[1]] + d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*sp[:x][k[1],j].out + d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    # objective needs to reflect the stage
    @stageobjective(sp, (1-pData.p0)^(n-1)*pData.p0*sp[:t][0].out);

    # add uncertainty
    SDDP.parameterize(sp, dDn, [item.prDis for item in dDn]) do ω
        JuMP.fix(H, ω.H)
        for i in pData.II
            if i != 0
                JuMP.fix(d[i],ω.d[i])
            end
        end
    end
end
