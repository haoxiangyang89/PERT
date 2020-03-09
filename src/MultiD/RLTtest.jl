using SDDP, Gurobi, JuMP;
using Distributions,HDF5,JLD,DelimitedFiles;
include("header.jl");

# initiate the data
filePath = "/Users/haoxiangyang/Desktop/PERT_tests/current/14";

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

    @constraint(mp, durationConstr1[k in pData.K, ω in Ω], t[k[2],ω] - t[k[1],ω] >= pData.D[k[1]] + disData[ω].d[k[1]]*G[k[1],ω]
                      - sum(pData.D[k[1]]*pData.eff[k[1]][j]*x[k[1],j,ω] + disData[ω].d[k[1]]*pData.eff[k[1]][j]*xG[k[1],j,ω] for j in pData.Ji[k[1]]));
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

    @constraint(mp, xGlinear1[i in pData.II, j in pData.Ji[i], ω in Ω], xG[i,j,ω] <= G[i,ω]);
    @constraint(mp, xGlinear2[i in pData.II, j in pData.Ji[i], ω in Ω], xG[i,j,ω] <= x[i,j,ω]);
    @constraint(mp, xGlinear3[i in pData.II, j in pData.Ji[i], ω in Ω], xG[i,j,ω] >= x[i,j,ω] - 1 + G[i,ω]);
    @constraint(mp, xGhatlinear1[i in pData.II, j in pData.Ji[i], ω in Ω], xGhat[i,j,ω] <= G[i,ω]);
    @constraint(mp, xGhatlinear2[i in pData.II, j in pData.Ji[i], ω in Ω], xGhat[i,j,ω] <= x0[i,j]);
    @constraint(mp, xGhatlinear3[i in pData.II, j in pData.Ji[i], ω in Ω], xGhat[i,j,ω] >= x0[i,j] - 1 + G[i,ω]);

    @objective(mp,Min,pData.p0*t0[0] + sum(disData[ω].prDis*t[0,ω] for ω in Ω));

    optimize!(mp);
    text = value.(mp[:t0]);
    xext = value.(mp[:x0]);
    gext = value.(mp[:G]);
    fext = objective_value(mp);
    return text,xext,fext,gext,mp;
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
tLRT,xLRT,fLRT,gLRT,mpLRT = extPERT_RLT(pData,disData,Ω);
