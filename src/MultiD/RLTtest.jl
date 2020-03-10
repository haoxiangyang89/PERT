using SDDP, Gurobi, JuMP;
using Distributions,HDF5,JLD,DelimitedFiles;
include("header.jl");

# initiate the data
#filePath = "/Users/haoxiangyang/Desktop/PERT_tests/current/14";
filePath = "/home/haoxiangyan/PERT_tests/current/14";   # for Darwin

# test RLT simple case function
function simpleRLTtest(that)
    sp = Model(Gurobi.Optimizer);
    @variable(sp, 0 <= y <= 1);
    @variable(sp, 0 <= G <= 1);

    @constraint(sp, y >= 0.5*G);
    @constraint(sp, y <= G);
    @constraint(sp, y <= that);
    @constraint(sp, y >= that + G - 1);
    @constraint(sp, that - y <= 0.5 - 0.5*G);
    @constraint(sp, that - y >= 0.5*G - 0.5);

    @objective(sp, Min, that + 0.5*(1 - G));
    optimize!(sp);
    return objective_value(sp);
end

# two-stage PERT extensive formulation
function extPERT(pData,disData,Ω,TL = Inf,relax = false)
    M = Dict();
    for ω in Ω
        M[ω] = sum(max(pData.D[i],pData.D[i]+disData[ω].d[i]) for i in pData.II if i != 0);
    end

    mp = Model(Gurobi.Optimizer);
    set_optimizer_attributes(mp, "IntFeasTol" => 1e-9, "TimeLimit" => TL);
    # mp = Model(solver = CplexSolver(CPX_PARAM_EPRHS = 1e-7,CPX_PARAM_EPINT = 1e-7,CPX_PARAM_TILIM = TL));
    @variable(mp,t0[i in pData.II] >= 0);
    @variable(mp,0 <= x0[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(mp,t[i in pData.II, ω in Ω] >= 0);
    @variable(mp,0 <= x[i in pData.II, j in pData.Ji[i], ω in Ω] <= 1);
    if relax
        @variable(mp,0 <= G[i in pData.II, ω in Ω] <= 1);
    else
        @variable(mp,G[i in pData.II, ω in Ω], Bin);
    end
    @variable(mp,0 <= s[i in pData.II, j in pData.Ji[i], ω in Ω] <= 1);

    @constraint(mp, FConstr[i in pData.II, ω in Ω], disData[ω].H - (1 - G[i,ω])*M[ω] <= t0[i]);
    @constraint(mp, GConstr[i in pData.II, ω in Ω], disData[ω].H + G[i,ω]*M[ω] >= t0[i]);
    @constraint(mp, tConstr1[i in pData.II, ω in Ω], t[i,ω] + G[i,ω]*M[ω] >= t0[i]);
    @constraint(mp, tConstr2[i in pData.II, ω in Ω], t[i,ω] - G[i,ω]*M[ω] <= t0[i]);
    @constraint(mp, tConstr3[i in pData.II, ω in Ω], t[i,ω] >= disData[ω].H * G[i,ω]);
    @constraint(mp, xConstr1[i in pData.II, j in pData.Ji[i], ω in Ω], x[i,j,ω] + G[i,ω] >= x0[i,j]);
    @constraint(mp, xConstr2[i in pData.II, j in pData.Ji[i], ω in Ω], x[i,j,ω] - G[i,ω] <= x0[i,j]);
    @constraint(mp, durationConstr1[k in pData.K, ω in Ω], t[k[2],ω] - t[k[1],ω] >= pData.D[k[1]] + disData[ω].d[k[1]]*G[k[1],ω]
                      - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j,ω] + disData[ω].d[k[1]]*pData.eff[k[1]][j]*s[k[1],j,ω] for j in pData.Ji[k[1]]));
    @constraint(mp, durationConstr2[k in pData.K], t0[k[2]] - t0[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x0[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II,ω in Ω], sum(x[i,j,ω] for j in pData.Ji[i]) <= 1);
    @constraint(mp, xConstr0[i in pData.II], sum(x0[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(mp, GGcons1[i in pData.II, ω in 1:(length(Ω) - 1)], G[i,ω] >= G[i,ω + 1]);
    @constraint(mp, GGcons2[k in pData.K, ω in Ω], G[k[2],ω] >= G[k[1],ω]);
    @constraint(mp, budgetConstr[ω in Ω], sum(sum(pData.b[i][j]*x[i,j,ω] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, budgetConstr0, sum(sum(pData.b[i][j]*x0[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, Slinear1[i in pData.II, j in pData.Ji[i], ω in Ω], s[i,j,ω] <= G[i,ω]);
    @constraint(mp, Slinear2[i in pData.II, j in pData.Ji[i], ω in Ω], s[i,j,ω] <= x[i,j,ω]);
    @constraint(mp, Slinear3[i in pData.II, j in pData.Ji[i], ω in Ω], s[i,j,ω] >= x[i,j,ω] - 1 + G[i,ω]);

    @objective(mp,Min,pData.p0*t0[0] + sum(disData[ω].prDis*t[0,ω] for ω in Ω));

    optimize!(mp);
    text = value.(mp[:t0]);
    xext = value.(mp[:x0]);
    text1 = value.(mp[:t]);
    xext1 = value.(mp[:x]);
    gext = value.(mp[:G]);
    fext = objective_value(mp);
    return text,xext,fext,text1,xext1,gext,mp;
end

# IP subs
function subIntC(pData,dDω,xhat,that,M = 999999,returnOpt = 0)
    # solve the MIP recourse problem
    sp = Model(Gurobi.Optimizer);
    set_optimizer_attributes(sp, "OutputFlag" => 0,"Threads" => 1);
    @variable(sp, 0 <= x[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, 0 <= s[i in pData.II,j in pData.Ji[i]] <= 1);
    @variable(sp, t[i in pData.II] >= 0);
    @variable(sp, G[i in pData.II], Bin);

    # add the basic sub problem constraints
    @constraint(sp, FCons[i in pData.II],dDω.H - (1 - G[i])*M <= that[i]);
    @constraint(sp, GCons[i in pData.II],dDω.H + G[i]*M >= that[i]);

    @constraint(sp, tGcons1[i in pData.II], t[i] <= that[i] + M*G[i]);
    @constraint(sp, tGcons2[i in pData.II], t[i] >= that[i] - M*G[i]);
    @constraint(sp, boundT[i in pData.II], t[i] >= dDω.H*G[i]);
    @constraint(sp, xGcons1[i in pData.II,j in pData.Ji[i]], x[i,j] <= xhat[i,j] + G[i]);
    @constraint(sp, xGcons2[i in pData.II,j in pData.Ji[i]], x[i,j] >= xhat[i,j] - G[i]);

    @constraint(sp, budgetConstr, sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(sp, xConstr[i in pData.II], sum(x[i,j] for j in pData.Ji[i]) <= 1);
    # linearize the bilinear term of x[i,j]*G[i]
    @constraint(sp, xGlin1[i in pData.II, j in pData.Ji[i]], s[i,j] <= G[i]);
    @constraint(sp, xGlin2[i in pData.II, j in pData.Ji[i]], s[i,j] <= x[i,j]);
    @constraint(sp, xGlin3[i in pData.II, j in pData.Ji[i]], s[i,j] >= x[i,j] - 1 + G[i]);

    @constraint(sp, durationConstr[k in pData.K], t[k[2]] - t[k[1]] >= pData.D[k[1]] + dDω.d[k[1]]*G[k[1]]
        - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] + dDω.d[k[1]]*pData.eff[k[1]][j]*s[k[1],j] for j in pData.Ji[k[1]]));

    @objective(sp, Min, t[0]);
    if returnOpt == 0
        optimize!(sp);
        return objective_value(sp);
    else
        return sp;
    end
end

# two-stage PERT RLT relaxation formulation
function extPERT_RLT(pData,disData,Ω,TL = Inf)
    M = Dict();
    for ω in Ω
        M[ω] = sum(max(pData.D[i],pData.D[i]+disData[ω].d[i]) for i in pData.II if i != 0);
    end

    mp = Model(Gurobi.Optimizer);
    set_optimizer_attributes(mp, "IntFeasTol" => 1e-6, "TimeLimit" => TL);
    # mp = Model(solver = CplexSolver(CPX_PARAM_EPRHS = 1e-7,CPX_PARAM_EPINT = 1e-7,CPX_PARAM_TILIM = TL));
    @variable(mp,t0[i in pData.II] >= 0);
    @variable(mp,0 <= x0[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(mp,t[i in pData.II, ω in Ω] >= 0);
    @variable(mp,0 <= x[i in pData.II, j in pData.Ji[i], ω in Ω] <= 1);
    @variable(mp,0 <= G[i in pData.II, ω in Ω] <= 1);
    @variable(mp,tG[i in pData.II, ω in Ω] >= 0);
    @variable(mp,0 <= xG[i in pData.II, j in pData.Ji[i], ω in Ω] <= 1);
    @variable(mp,tGhat[i in pData.II, ω in Ω] >= 0);
    @variable(mp,0 <= xGhat[i in pData.II, j in pData.Ji[i], ω in Ω] <= 1);
    @variable(mp,tGR1[k in pData.K, ω in Ω] >= 0);
    @variable(mp,tGR2[k in pData.K, ω in Ω] >= 0);
    @variable(mp,0 <= GG[k in pData.K, ω in Ω] <= 1);
    @variable(mp,0 <= xGR[k in pData.K, j in pData.Ji[k[1]], ω in Ω] <= 1);
    @variable(mp,0 <= xGG[k in pData.K, j in pData.Ji[k[1]], ω in Ω] <= 1);

    @constraint(mp, FConstr1[i in pData.II, ω in Ω], disData[ω].H*G[i,ω] <= tGhat[i,ω]);
    @constraint(mp, FConstr2[i in pData.II, ω in Ω], (disData[ω].H - M[ω])*(1 - G[i,ω]) <= t0[i] - tGhat[i,ω]);
    @constraint(mp, GConstr1[i in pData.II, ω in Ω], disData[ω].H*G[i,ω] + G[i,ω]*M[ω] >= tGhat[i,ω]);
    @constraint(mp, GConstr2[i in pData.II, ω in Ω], disData[ω].H*(1 - G[i,ω]) >= t0[i] - tGhat[i,ω]);

    @constraint(mp, tConstr11[i in pData.II, ω in Ω], tG[i,ω] + G[i,ω]*M[ω] >= tGhat[i,ω]);
    @constraint(mp, tConstr12[i in pData.II, ω in Ω], t[i,ω] - tG[i,ω] == t0[i] - tGhat[i,ω]);

    @constraint(mp, tConstr21[i in pData.II, ω in Ω], tG[i,ω] - G[i,ω]*M[ω] <= tGhat[i,ω]);

    @constraint(mp, tConstr31[i in pData.II, ω in Ω], tG[i,ω] >= disData[ω].H * G[i,ω]);
    @constraint(mp, tConstr32[i in pData.II, ω in Ω], t[i,ω] >= disData[ω].H * G[i,ω]);

    @constraint(mp, xConstr11[i in pData.II, j in pData.Ji[i], ω in Ω], xG[i,j,ω] + G[i,ω] >= xGhat[i,j,ω]);
    @constraint(mp, xConstr12[i in pData.II, j in pData.Ji[i], ω in Ω], x[i,j,ω] - xG[i,j,ω] == x0[i,j] - xGhat[i,j,ω]);

    @constraint(mp, xConstr21[i in pData.II, j in pData.Ji[i], ω in Ω], xG[i,j,ω] - G[i,ω] <= xGhat[i,j,ω]);

    @constraint(mp, durationConstr11[k in pData.K, ω in Ω], tGR1[k,ω] - tG[k[1],ω] >= pData.D[k[1]]*G[k[1],ω] + disData[ω].d[k[1]]*G[k[1],ω]
                      - sum(pData.D[k[1]]*pData.eff[k[1]][j]*xG[k[1],j,ω] + disData[ω].d[k[1]]*pData.eff[k[1]][j]*xG[k[1],j,ω] for j in pData.Ji[k[1]]));
    @constraint(mp, durationConstr12[k in pData.K, ω in Ω], t[k[2],ω] - tGR1[k,ω] + tG[k[1],ω] - t[k[1],ω] >= pData.D[k[1]]*(1 - G[k[1],ω])
                      - sum(pData.D[k[1]]*pData.eff[k[1]][j]*(x[k[1],j,ω] - xG[k[1],j,ω]) for j in pData.Ji[k[1]]));
    @constraint(mp, durationConstr13[k in pData.K, ω in Ω], tG[k[2],ω] - tGR2[k,ω] >= pData.D[k[1]]*G[k[2],ω] + disData[ω].d[k[1]]*GG[k,ω]
                      - sum(pData.D[k[1]]*pData.eff[k[1]][j]*xGR[k,j,ω] + disData[ω].d[k[1]]*pData.eff[k[1]][j]*xGG[k,j,ω] for j in pData.Ji[k[1]]));
    @constraint(mp, durationConstr14[k in pData.K, ω in Ω], t[k[2],ω] - tG[k[2],ω] + tGR2[k,ω] - t[k[1],ω] >= pData.D[k[1]]*(1 - G[k[2],ω]) + disData[ω].d[k[1]]*(G[k[1],ω] - GG[k,ω])
                      - sum(pData.D[k[1]]*pData.eff[k[1]][j]*(x[k[1],j,ω] - xGR[k,j,ω]) + disData[ω].d[k[1]]*pData.eff[k[1]][j]*(xG[k[1],j,ω] - xGG[k,j,ω]) for j in pData.Ji[k[1]]));

    @constraint(mp, durationConstr2[k in pData.K], t0[k[2]] - t0[k[1]] >= pData.D[k[1]]*(1 - sum(pData.eff[k[1]][j]*x0[k[1],j] for j in pData.Ji[k[1]])));
    @constraint(mp, xConstr[i in pData.II,ω in Ω], sum(x[i,j,ω] for j in pData.Ji[i]) <= 1);
    @constraint(mp, xConstr0[i in pData.II], sum(x0[i,j] for j in pData.Ji[i]) <= 1);
    @constraint(mp, GGcons1[i in pData.II, ω in 1:(length(Ω) - 1)], G[i,ω] >= G[i,ω + 1]);
    @constraint(mp, GGcons2[k in pData.K, ω in Ω], G[k[2],ω] >= G[k[1],ω]);
    @constraint(mp, budgetConstr[ω in Ω], sum(sum(pData.b[i][j]*x[i,j,ω] for j in pData.Ji[i]) for i in pData.II) <= pData.B);
    @constraint(mp, budgetConstr0, sum(sum(pData.b[i][j]*x0[i,j] for j in pData.Ji[i]) for i in pData.II) <= pData.B);

    @constraint(mp, tGlinear1[i in pData.II, ω in Ω], tG[i,ω] <= M[ω]*G[i,ω]);
    @constraint(mp, tGlinear2[i in pData.II, ω in Ω], tG[i,ω] <= t[i,ω]);
    @constraint(mp, tGlinear3[i in pData.II, ω in Ω], tG[i,ω] >= t[i,ω] + M[ω]*(G[i,ω] - 1));
    @constraint(mp, tGhatlinear1[i in pData.II, ω in Ω], tGhat[i,ω] <= M[ω]*G[i,ω]);
    @constraint(mp, tGhatlinear2[i in pData.II, ω in Ω], tGhat[i,ω] <= t0[i]);
    @constraint(mp, tGhatlinear3[i in pData.II, ω in Ω], tGhat[i,ω] >= t0[i] + M[ω]*(G[i,ω] - 1));

    @constraint(mp, tGRlinear11[k in pData.K, ω in Ω], tGR1[k,ω] <= M[ω]*G[k[1],ω]);
    @constraint(mp, tGRlinear12[k in pData.K, ω in Ω], tGR1[k,ω] <= t[k[2],ω]);
    @constraint(mp, tGRlinear13[k in pData.K, ω in Ω], tGR1[k,ω] >= t[k[2],ω] + M[ω]*(G[k[1],ω] - 1));
    @constraint(mp, tGRlinear21[k in pData.K, ω in Ω], tGR2[k,ω] <= M[ω]*G[k[2],ω]);
    @constraint(mp, tGRlinear22[k in pData.K, ω in Ω], tGR2[k,ω] <= t[k[1],ω]);
    @constraint(mp, tGRlinear23[k in pData.K, ω in Ω], tGR2[k,ω] >= t[k[1],ω] + M[ω]*(G[k[2],ω] - 1));

    @constraint(mp, GGlinear1[k in pData.K, ω in Ω], GG[k,ω] <= G[k[1],ω]);
    @constraint(mp, GGlinear2[k in pData.K, ω in Ω], GG[k,ω] <= G[k[2],ω]);
    @constraint(mp, GGlinear3[k in pData.K, ω in Ω], GG[k,ω] >= G[k[1],ω] + G[k[2],ω] - 1);

    @constraint(mp, xGlinear1[i in pData.II, j in pData.Ji[i], ω in Ω], xG[i,j,ω] <= G[i,ω]);
    @constraint(mp, xGlinear2[i in pData.II, j in pData.Ji[i], ω in Ω], xG[i,j,ω] <= x[i,j,ω]);
    @constraint(mp, xGlinear3[i in pData.II, j in pData.Ji[i], ω in Ω], xG[i,j,ω] >= x[i,j,ω] - 1 + G[i,ω]);
    @constraint(mp, xGhatlinear1[i in pData.II, j in pData.Ji[i], ω in Ω], xGhat[i,j,ω] <= G[i,ω]);
    @constraint(mp, xGhatlinear2[i in pData.II, j in pData.Ji[i], ω in Ω], xGhat[i,j,ω] <= x0[i,j]);
    @constraint(mp, xGhatlinear3[i in pData.II, j in pData.Ji[i], ω in Ω], xGhat[i,j,ω] >= x0[i,j] - 1 + G[i,ω]);

    @constraint(mp, xGRlinear1[k in pData.K, j in pData.Ji[k[1]], ω in Ω], xGR[k,j,ω] <= G[k[2],ω]);
    @constraint(mp, xGRlinear2[k in pData.K, j in pData.Ji[k[1]], ω in Ω], xGR[k,j,ω] <= x[k[1],j,ω]);
    @constraint(mp, xGRlinear3[k in pData.K, j in pData.Ji[k[1]], ω in Ω], xGR[k,j,ω] >= x[k[1],j,ω] - 1 + G[k[2],ω]);
    @constraint(mp, xGGlinear1[k in pData.K, j in pData.Ji[k[1]], ω in Ω], xGG[k,j,ω] <= GG[k,ω]);
    @constraint(mp, xGGlinear2[k in pData.K, j in pData.Ji[k[1]], ω in Ω], xGG[k,j,ω] <= x[k[1],j,ω]);
    @constraint(mp, xGGlinear3[k in pData.K, j in pData.Ji[k[1]], ω in Ω], xGG[k,j,ω] >= x[k[1],j,ω] - 1 + GG[k,ω]);

    @objective(mp,Min,pData.p0*t0[0] + sum(disData[ω].prDis*t[0,ω] for ω in Ω));

    optimize!(mp);
    text = value.(mp[:t0]);
    xext = value.(mp[:x0]);
    gext = value.(mp[:G]);
    fext = objective_value(mp);
    return text,xext,fext,gext,mp;
end

# solve the lagrangian relaxation on top of the RLT
function RLT_LR_sub(pData,dDω,λDict,μ,t0,x0,Mω)
    sp = Model(Gurobi.Optimizer);
    set_optimizer_attributes(sp, "IntFeasTol" => 1e-6, "OutputFlag" => 0);
    @variable(sp, t[i in pData.II] >= 0);
    @variable(sp, 0 <= x[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(sp, 0 <= G[i in pData.II] <= 1);

    @variable(sp, tG[i in pData.II] >= 0);
    @variable(sp, 0 <= xG[i in pData.II, j in pData.Ji[i]] <= 1);
    @variable(sp, tGhat[i in pData.II] >= 0);
    @variable(sp, 0 <= xGhat[i in pData.II, j in pData.Ji[i]] <= 1);

    @constraint(sp, FConstr1[i in pData.II], dDω.H*G[i] <= tGhat[i]);
    @constraint(sp, FConstr2[i in pData.II], (dDω.H - Mω)*(1 - G[i]) <= t0[i] - tGhat[i]);
    @constraint(sp, GConstr1[i in pData.II], dDω.H*G[i] + G[i]*Mω >= tGhat[i]);
    @constraint(sp, GConstr2[i in pData.II], dDω.H*(1 - G[i]) >= t0[i] - tGhat[i]);
    @constraint(sp, tConstr11[i in pData.II], tG[i] + G[i]*Mω >= tGhat[i]);
    @constraint(sp, tConstr12[i in pData.II], t[i] - tG[i] == t0[i] - tGhat[i]);
    @constraint(sp, tConstr21[i in pData.II], tG[i] - G[i]*Mω <= tGhat[i]);
    @constraint(sp, tConstr31[i in pData.II], tG[i] >= dDω.H * G[i]);
    @constraint(sp, tConstr32[i in pData.II], t[i] >= dDω.H * G[i]);
    @constraint(sp, xConstr11[i in pData.II, j in pData.Ji[i]], xG[i,j] + G[i] >= xGhat[i,j]);
    @constraint(sp, xConstr12[i in pData.II, j in pData.Ji[i]], x[i,j] - xG[i,j] == x0[i,j] - xGhat[i,j]);
    @constraint(sp, xConstr21[i in pData.II, j in pData.Ji[i]], xG[i,j] - G[i] <= xGhat[i,j]);

    @constraint(sp, tGlinear1[i in pData.II], tG[i] <= Mω*G[i]);
    @constraint(sp, tGlinear2[i in pData.II], tG[i] <= t[i]);
    @constraint(sp, tGlinear3[i in pData.II], tG[i] >= t[i] + Mω*(G[i] - 1));
    @constraint(sp, tGhatlinear1[i in pData.II], tGhat[i] <= Mω*G[i]);
    @constraint(sp, tGhatlinear2[i in pData.II], tGhat[i] <= t0[i]);
    @constraint(sp, tGhatlinear3[i in pData.II], tGhat[i] >= t0[i] + Mω*(G[i] - 1));

    @constraint(sp, xGlinear1[i in pData.II, j in pData.Ji[i]], xG[i,j] <= G[i]);
    @constraint(sp, xGlinear2[i in pData.II, j in pData.Ji[i]], xG[i,j] <= x[i,j]);
    @constraint(sp, xGlinear3[i in pData.II, j in pData.Ji[i]], xG[i,j] >= x[i,j] - 1 + G[i]);
    @constraint(sp, xGhatlinear1[i in pData.II, j in pData.Ji[i]], xGhat[i,j] <= G[i]);
    @constraint(sp, xGhatlinear2[i in pData.II, j in pData.Ji[i]], xGhat[i,j] <= x0[i,j]);
    @constraint(sp, xGhatlinear3[i in pData.II, j in pData.Ji[i]], xGhat[i,j] >= x0[i,j] - 1 + G[i]);
    
    # including the Lagrangian relaxation term
    @expression(sp, λCoeff[k in pData.K], pData.D[k[1]] + dDω.d[k[1]]*G[k[1]] - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j] +
        dDω.d[k[1]]*pData.eff[k[1]][j]*xG[k[1],j] for j in pData.Ji[k[1]]) + t[k[1]] - t[k[2]]);
    @expression(sp, μCoeff, -pData.B + sum(sum(pData.b[i][j]*x[i,j] for j in pData.Ji[i]) for i in pData.II));
    @objective(sp, Min, t[0] + sum(λDict[k]*λCoeff[k] for k in pData.K) + μ*μCoeff);

    optimize!(sp);
    tLD = value.(sp[:t]);
    xLD = value.(sp[:x]);
    gLD = value.(sp[:G]);
    λCoeffLD = value.(sp[:λCoeff]);
    μCoeffLD = value.(sp[:μCoeff]);
    fLD = objective_value(sp);
    return tLD,xLD,fLD,gLD,λCoeffLD,μCoeffLD;
end

# Lagrangian relaxation process
function RLT_LR_Proc(pData,disData,Ω,λ0,μ0)
    M = Dict();
    for ω in Ω
        M[ω] = sum(max(pData.D[i],pData.D[i]+disData[ω].d[i]) for i in pData.II if i != 0);
    end

    # solve the master, obtain a t0 and x0

    for ω in Ω
        # initial state of LD multiplier
        λDict = deepcopy(λ0[ω]);
        μ = μ0[ω];
        global bestf = -Inf;
        global bestResult = [];
        for N in 1:1000
            # solve the Lagrangian dual
            tLD,xLD,fLD,gLD,λCoeffLD,μCoeffLD = RLT_LR_sub(pData,disData[ω],λDict,μ,t0,x0,M[ω]);
            println("============================= Iteration $(N), Current LB $(fLD), Best LB $(bestf) =============================");
            if fLD > bestf
                global bestf = fLD;
                global bestResult = [tLD,xLD,fLD,gLD,deepcopy(λDict),μ];
            end
            # update the multiplier
            for k in pData.K
                if abs(λCoeffLD[k]) > 1e-8
                    if λDict[k] + λCoeffLD[k]*α < 0
                        global λDict[k] = 0;
                    else
                        global λDict[k] += λCoeffLD[k]*α;
                    end
                    # global λDict[k] += λCoeffLD[k]*α;
                end
            end
            if abs(μCoeffLD) > 1e-8
                if μ + μCoeffLD*α < 0
                    global μ = 0;
                else
                    global μ += μCoeffLD*α;
                end
                # global μ += μCoeffLD*α;
            end
        end
    end
    # generate cuts from sub and feed back to the master
end


# ==============================================================================
# a simple RLT test for one indicator, one t setting, tight approximation of the value function
tList = 0:0.01:1;
objList = [];
for t in tList
    push!(objList,simpleRLTtest(t));
end

# generate two stage data
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,1,999999);
disData,Ω = autoUGen(nameH, Hparams, nameD, dparams, 100, 1);
disData = orderdisData(disData,Ω);
# IP extensive form and LP relaxation
text,xext,fext,text1,xext1,gext,mp = extPERT(pData,disData,Ω);
tLPR,xLPR,fLPR,gLPR,mpLPR = extPERT(pData,disData,Ω,Inf,true);
tRLT,xRLT,fRLT,gRLT,mpRLT = extPERT_RLT(pData,disData,Ω);

# LD method
λ0 = Dict();
μ0 = Dict();

for ω in Ω
    λ0[ω] = Dict();
    for k in pData.K
        λ0[ω][k] = 0;
    end
    μ0[ω] = 0;
end
t0 = Dict();
x0 = Dict();
for i in pData.II
    t0[i] = text[i];
    for j in pData.Ji[i]
        x0[i,j] = xext[i,j];
    end
end
