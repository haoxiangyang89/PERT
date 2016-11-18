# This is a collections of all functions
using JuMP,Distributions,CPLEX;

function readIn(InputAdd)
    # no of activities
    data = readdlm(InputAdd,',');
    lI = data[1,1];
    II = 1:lI;
    lJ = data[1,1];
    J = 1:data[1,2];
    row = 3;

    # the precedence relationship
    G = [];
    while data[row,1]!=""
        push!(G,data[row,1:2]);
        row += 1;
    end
    row += 1;

    # the nominal duration/after disruption ratio information
    D = Dict();
    while data[row,1]!=""
        D[data[row,1]] = data[row,2];
        row += 1;
    end
    row += 1;

    # the budget and the cost to crash
    B = data[row,1];
    ee = Dict();
    for j in J
        ee[j] = data[row,Int64(j+1)];
    end
    row += 1;
    b = Dict();
    while data[row,1]!=""
        for j in J
            b[data[row,1],j] = data[row,Int64(j+1)];
        end
        row += 1;
    end

    # sample according to the distribution
    r = Dict();
    H = Dict();
    # no disruption in first scenario

    # first read in the distribution names and then the parameters
    row += 1;
    Hname = data[row,1];
    if Hname == "Exponential"
        μ = data[row,2];
        distrH = Exponential(μ);
    elseif Hname == "LogNormal"
        μ = data[row,2];
        σ = data[row,3];
        distrH = Normal(μ,σ);
    elseif Hname == "Gamma"
        α = data[row,2];
        θ = data[row,3];
        distrH = Gamma(α,θ);
    elseif Hname == "Uniform"
        la = data[row,2];
        ub = data[row,3];
        distrH = Uniform(la,ub+1e-10);
    end
    row += 1;
    rname = data[row,1];
    row += 1;
    if rname == "Uniform"
        la = Dict();
        ub = Dict();
        while data[row,1]!=""
            la[data[row,1]] = data[row,2];
            ub[data[row,1]] = data[row,3];
            row += 1;
        end
        disRList = Dict();
        for i in II
            r[i,1] = 1;
            disRList[i] = Uniform(la[i],ub[i]+1e-10);
        end
    elseif rname == "LogNormal"
        μ = Dict();
        σ = Dict();
        while data[row,1]!=""
            μ[data[row,1]] = data[row,2];
            σ[data[row,1]] = data[row,3];
            row += 1;
        end
        disRList = Dict();
        for i in II
            r[i,1] = 1;
            disRList[i] = LogNormal(μ,σ);
        end
    end

    # the number of scenarios
    row += 1;
    S = data[row,1] + 1;
    p = data[row,2];
    M = Dict();

    for s in 2:S
        H[s] = round(rand(distrH),4);
        while H[s] > 97
            H[s] = round(rand(distrH),4);
        end
        for i in II
            r[i,s] = round(rand(disRList[i]),4);
        end
        M[s] = sum(values(D))*maximum([r[i,s] for i in II]);
    end
    H[1] = maximum(values(M));
    M[1] = 2*sum(values(D));
    SS = 1:S;

    return D,r,H,b,B,ee,II,J,M,SS,G,distrH,disRList,p
end

function makeMaster(D,r,H,b,B,ee,II,JJ,SS,GG,dH,dR,p)
    # this function builds up the master program
    #env = Gurobi.Env();
    #setparams!(env;MIPGap = 1e-8,MIPFocus = 2);
    #solver = GurobiSolver(env);
    solver = CplexSolver(CPX_PARAM_THREADS = 8);

    # construct a master problem
    m = Model();
    @variables(m,begin
        θ[s in SS;s>=2] >= 0
        t[i in II] >= 0
        tN >= 0
        x[i in II, j in JJ], Bin
    end);

    @objective(m,Min,p*tN + (1-p)/(length(SS)-1)*sum{θ[s],s in SS;s>=2});
    @constraint(m,sepN[i in II],tN>=t[i]);
    @constraint(m,sepi[g in 1:length(GG)],t[GG[g][2]]-t[GG[g][1]]>=D[GG[g][1]]*(1-sum{ee[j]*x[GG[g][1],j],j in JJ}));
    @constraint(m,oneAction[i in II],sum{x[i,j],j in JJ} <= 1);
    @constraint(m,budget,sum{x[i,j]*b[i,j],i in II, j in JJ} <= B);

    return m
end

function obtainIs(xmaster,tmaster,D,rsub,Hsub,b,B,ee,II,JJ,GG)
    # obtain three categories of positions of events for the second stage
    I1 = [];
    I2 = [];
    I3 = [];
    tau = Dict{Any,Any}();
    for i in II
        if tmaster[i] >= Hsub
            push!(I3,i);
            tau[i] = D[i]*rsub[i];
        else
            if tmaster[i] + (1 - sum([ee[j]*xmaster[i,j] for j in JJ]))*D[i] > Hsub
                push!(I2,i);
                tau[i] = (Hsub - tmaster[i]) + rsub[i]*(tmaster[i] + (1 - sum([ee[j]*xmaster[i,j] for j in JJ]))*D[i] - Hsub);
            else
                push!(I1,i);
                tau[i] = D[i]*(1-sum([ee[j]*xmaster[i,j] for j in JJ]));
            end
        end
    end
    return I1,I2,I3,tau
end

function subInt(xmaster,tmaster,D,tau,b,B,ee,II,I1,I2,I3,JJ,GG)
    # build the subproblem to obtain the upper bound
    msub = Model();
    @variables(msub,begin
        0 <= x1[i in I1, j in JJ] <= 1
        0 <= x2[i in I2, j in JJ] <= 1
        x3[i in I3, j in JJ], Bin
        te[i in II] >= 0
        tN >= 0
    end);
    @constraint(msub,bindingX1[i in I1, j in JJ],x1[i,j] == xmaster[i,j]);
    @constraint(msub,bindingX2[i in I2, j in JJ],x2[i,j] == xmaster[i,j]);
    @constraint(msub,binding1[i in I1],te[i] == tmaster[i]);
    @constraint(msub,binding2[i in I2],te[i] == tmaster[i]);
    @constraint(msub,sepN[i in II], tN >= te[i]);
    @constraint(msub,sepi1[g in 1:length(GG);GG[g][1] in I1], te[GG[g][2]] - te[GG[g][1]] >= tau[GG[g][1]]);
    @constraint(msub,sepi2[g in 1:length(GG);GG[g][1] in I2], te[GG[g][2]] - te[GG[g][1]] >= tau[GG[g][1]]);
    @constraint(msub,sepi3[g in 1:length(GG);GG[g][1] in I3], te[GG[g][2]] - te[GG[g][1]] >= tau[GG[g][1]]*(1-sum{ee[j]*x3[GG[g][1],j], j in JJ}));
    @constraint(msub,oneAction1[i in I1], sum{x1[i,j], j in JJ} <= 1);
    @constraint(msub,oneAction2[i in I2], sum{x2[i,j], j in JJ} <= 1);
    @constraint(msub,oneAction3[i in I3], sum{x3[i,j], j in JJ} <= 1);
    @constraint(msub,budget,sum{x1[i,j]*b[i,j], i in I1, j in JJ} +
      sum{x2[i,j]*b[i,j], i in I2, j in JJ} + sum{x3[i,j]*b[i,j], i in I3, j in JJ}<= B);
    @objective(msub,Min,tN);

    return msub;
end

function subLagP(xmaster,tmaster,D,Hsub,tau,b,B,ee,II,I1,I2,I3,JJ,GG,lambda1,lambda2,ppi)
    # build the subproblem for easy model of lagrangian relaxation
    mlagP = Model();
    @variables(mlagP,begin
        0 <= xr1[i in I1, j in JJ] <= 1
        0 <= xr2[i in I2, j in JJ] <= 1
        xr3[i in I3, j in JJ], Bin
        ter[i in II] >= 0
        0 <= tNr <= sum(values(tau));
    end);
    @constraint(mlagP,sepN[i in II], tNr >= ter[i]);
    @constraint(mlagP,sepi1[g in 1:length(GG);GG[g][1] in I1], ter[GG[g][2]] - ter[GG[g][1]] >= tau[GG[g][1]]);
    @constraint(mlagP,sepi2[g in 1:length(GG);GG[g][1] in I2], ter[GG[g][2]] - ter[GG[g][1]] >= tau[GG[g][1]]);
    @constraint(mlagP,sepi3[g in 1:length(GG);GG[g][1] in I3], ter[GG[g][2]] - ter[GG[g][1]] >= tau[GG[g][1]]*(1-sum{ee[j]*xr3[GG[g][1],j], j in JJ}));
    @constraint(mlagP,oneAction1[i in I1], sum{xr1[i,j], j in JJ} <= 1);
    @constraint(mlagP,oneAction2[i in I2], sum{xr2[i,j], j in JJ} <= 1);
    @constraint(mlagP,oneAction3[i in I3], sum{xr3[i,j], j in JJ} <= 1);
    @constraint(mlagP,budget,sum{xr1[i,j]*b[i,j], i in I1, j in JJ} +
      sum{xr2[i,j]*b[i,j], i in I2, j in JJ} + sum{xr3[i,j]*b[i,j], i in I3, j in JJ}<= B);
    @objective(mlagP,Min,tNr+sum{lambda1[i,j]*(xmaster[i,j] - xr1[i,j]),i in I1,j in JJ}
        +sum{lambda2[i,j]*(xmaster[i,j]-xr2[i,j]),i in I2,j in JJ}
        +sum{ppi[i]*(tmaster[i]-ter[i]),i in I1}+sum{ppi[i]*(tmaster[i]-ter[i]),i in I2});

    return mlagP;
end

function subLag(xmaster,tmaster,D,rsub,Hsub,M,b,B,ee,II,I1,I2,I3,JJ,GG,lambda1,lambda2,ppi)
    # build the subproblem for hard model of lagrangian relaxation
    mlag = Model();
    @variables(mlag,begin
        0 <= xre1[i in I1, j in JJ] <= 1
        0 <= xre2[i in I2, j in JJ] <= 1
        xre3[i in I3, j in JJ], Bin
        tere[i in II] >= 0
        0 <= tNre <= M
        u[i in II] >= 0
        z[i in II, j in JJ], Bin
        F[i in II], Bin
        F[i in II], Bin
        Ft[i in II], Bin
        G[i in II], Bin
        Gt[i in II], Bin
        S1[i in II], Bin
        S2[i in II], Bin
        S3[i in II], Bin
        Z1[i in II, j in JJ], Bin
        Z2[i in II, j in JJ], Bin
        Z3[i in II, j in JJ], Bin
        Y[i in II] >= 0
    end);
    @objective(mlag,Min,tNre+sum{lambda1[i,j]*(xmaster[i,j] - xre1[i,j]),i in I1,j in JJ}
        +sum{lambda2[i,j]*(xmaster[i,j]-xre2[i,j]),i in I2,j in JJ}
        +sum{ppi[i]*(tmaster[i]-tere[i]),i in I1}+sum{ppi[i]*(tmaster[i]-tere[i]),i in I2});
    @constraint(mlag,timeInd1[i in II],Hsub - F[i]*M <= tere[i]);
    @constraint(mlag,timeInd2[i in II],Hsub + G[i]*M >= tere[i]);
    @constraint(mlag,timeIndC[i in II],F[i] + G[i] == 1);
    @constraint(mlag,EndInd11[i in I1],Hsub - Ft[i]*M <= tere[i] + D[i]*(1 - sum{ee[j]*xre1[i,j],j in JJ}));
    @constraint(mlag,EndInd12[i in I2],Hsub - Ft[i]*M <= tere[i] + D[i]*(1 - sum{ee[j]*xre2[i,j],j in JJ}));
    @constraint(mlag,EndInd13[i in I3],Hsub - Ft[i]*M <= tere[i] + D[i]*(1 - sum{ee[j]*xre3[i,j],j in JJ}));
    @constraint(mlag,EndInd21[i in I1],Hsub + Gt[i]*M >= tere[i] + D[i]*(1 - sum{ee[j]*xre1[i,j],j in JJ}));
    @constraint(mlag,EndInd22[i in I2],Hsub + Gt[i]*M >= tere[i] + D[i]*(1 - sum{ee[j]*xre2[i,j],j in JJ}));
    @constraint(mlag,EndInd23[i in I3],Hsub + Gt[i]*M >= tere[i] + D[i]*(1 - sum{ee[j]*xre3[i,j],j in JJ}));
    @constraint(mlag,EndIndC[i in II],Ft[i] + Gt[i] == 1);
    @constraint(mlag,nonAT1[i in II],u[i] >= tere[i] - (1 - F[i])*M);
    @constraint(mlag,nonAT2[i in II],u[i] <= tere[i] + (1 - F[i])*M);
    @constraint(mlag,nonAX11[i in I1, j in JJ],z[i,j] >= xre1[i,j] - (1 - F[i]));
    @constraint(mlag,nonAX21[i in I1, j in JJ],z[i,j] <= xre1[i,j] + (1 - F[i]));
    @constraint(mlag,nonAX12[i in I2, j in JJ],z[i,j] >= xre2[i,j] - (1 - F[i]));
    @constraint(mlag,nonAX22[i in I2, j in JJ],z[i,j] <= xre2[i,j] + (1 - F[i]));
    @constraint(mlag,nonAX13[i in I3, j in JJ],z[i,j] >= xre3[i,j] - (1 - F[i]));
    @constraint(mlag,nonAX23[i in I3, j in JJ],z[i,j] <= xre3[i,j] + (1 - F[i]));
    @constraint(mlag,sepN[i in II],tNre >= u[i]);
    @constraint(mlag,sepi[g in 1:length(GG)],u[GG[g][2]] - u[GG[g][1]] >=
    (rsub[GG[g][1]] - 1)*Y[GG[g][1]] + (rsub[GG[g][1]]*D[GG[g][1]]+Hsub*(1-rsub[GG[g][1]]))*S1[GG[g][1]]
    - (rsub[GG[g][1]]*D[GG[g][1]]) * sum{(ee[j]*Z1[GG[g][1],j]),j in JJ}
    + D[GG[g][1]]*(S2[GG[g][1]] - sum{ee[j]*Z2[GG[g][1],j],j in JJ})
    + D[GG[g][1]]*rsub[GG[g][1]]*(S3[GG[g][1]] - sum{ee[j]*Z3[GG[g][1],j],j in JJ})
    );
    @constraint(mlag,oneAction[i in II],sum{z[i,j],j in JJ} <= 1);
    @constraint(mlag,budget,sum{z[i,j]*b[i,j] ,i in II, j in JJ} <= B);
    @constraint(mlag,timecri1[i in II], F[i] + Gt[i] <= 1 + S1[i]);
    @constraint(mlag,timecri2[i in II], F[i] + Ft[i] <= 1 + S2[i]);
    @constraint(mlag,timecri3[i in II], G[i] + Gt[i] <= 1 + S3[i]);
    @constraint(mlag,timecriBind[i in II], S1[i] + S2[i] + S3[i] == 1);
    @constraint(mlag,linearZ11[i in II, j in JJ], Z1[i,j] <= S1[i]);
    @constraint(mlag,linearZ12[i in II, j in JJ], Z1[i,j] <= z[i,j]);
    @constraint(mlag,linearZ13[i in II, j in JJ], Z1[i,j] >= S1[i] + z[i,j] - 1);
    @constraint(mlag,linearZ21[i in II, j in JJ], Z2[i,j] <= S2[i]);
    @constraint(mlag,linearZ22[i in II, j in JJ], Z2[i,j] <= z[i,j]);
    @constraint(mlag,linearZ23[i in II, j in JJ], Z2[i,j] >= S2[i] + z[i,j] - 1);
    @constraint(mlag,linearZ31[i in II, j in JJ], Z3[i,j] <= S3[i]);
    @constraint(mlag,linearZ32[i in II, j in JJ], Z3[i,j] <= z[i,j]);
    @constraint(mlag,linearZ33[i in II, j in JJ], Z3[i,j] >= S3[i] + z[i,j] - 1);
    @constraint(mlag,linearY1[i in II], Y[i] <= M*S1[i]);
    @constraint(mlag,linearY2[i in II], Y[i] <= u[i]);
    @constraint(mlag,linearY3[i in II], Y[i] >= u[i] - M*(1-S1[i]));

    return mlag;
end

function fullExt(D,r,H,b,B,ee,II,JJ,SS,GG,dH,dR,p,M)
  m = Model();
  # set up the solver as Cplex
  solver = CplexSolver(CPX_PARAM_THREADS = 8);

  @variables(m,begin
      t[i in II, s in SS] >= 0
      tN[s in SS] >= 0
      F[i in II, s in SS], Bin
      Ft[i in II, s in SS], Bin
      G[i in II, s in SS], Bin
      Gt[i in II, s in SS], Bin
      x[i in II, j in JJ, s in SS], Bin
      S1[i in II, s in SS], Bin
      S2[i in II, s in SS], Bin
      S3[i in II, s in SS], Bin
      Z1[i in II, j in JJ, s in SS], Bin
      Z2[i in II, j in JJ, s in SS], Bin
      Z3[i in II, j in JJ, s in SS], Bin
      Y[i in II, s in SS] >= 0
  end)

  @objective(m,Min,(1-p)/(length(SS)-1)*sum{tN[s], s in SS;s!=1}+p*(tN[1]));
  @constraint(m,sepi1[g in 1:length(GG)],t[GG[g][2],1] - t[GG[g][1],1] >= D[GG[g][1]]*(1 - sum{ee[j]*x[GG[g][1],j,1],j in JJ}));
  @constraint(m,sepi[g in 1:length(GG),s in SS;s>=2],
    t[GG[g][2],s] - t[GG[g][1],s] >=
    (r[GG[g][1],s] - 1)*Y[GG[g][1],s] + (r[GG[g][1],s]*D[GG[g][1]]+H[s]*(1-r[GG[g][1],s]))*S1[GG[g][1],s]
    - (r[GG[g][1],s]*D[GG[g][1]]) * sum{(ee[j]*Z1[GG[g][1],j,s]),j in JJ}
    + D[GG[g][1]]*(S2[GG[g][1],s] - sum{ee[j]*Z2[GG[g][1],j,s],j in JJ})
    + D[GG[g][1]]*r[GG[g][1],s]*(S3[GG[g][1],s] - sum{ee[j]*Z3[GG[g][1],j,s],j in JJ})
    );
  @constraint(m,sepN[i in II, s in SS],tN[s] >= t[i,s]);
  @constraint(m,timeInd1[i in II, s in SS;s>=2],t[i,1] >= H[s]-F[i,s]*M[s]);
  @constraint(m,timeInd2[i in II, s in SS;s>=2],t[i,1] <= H[s]+G[i,s]*M[s]);
  @constraint(m,timeIndC[i in II, s in SS;s>=2],G[i,s] + F[i,s] == 1);
  @constraint(m,EndInd1[i in II, s in SS;s>=2],t[i,1] + D[i]*(1-sum{ee[j]*x[i,j,s],j in JJ}) >= H[s]-Ft[i,s]*M[s]);
  @constraint(m,EndInd2[i in II, s in SS;s>=2],t[i,1] + D[i]*(1-sum{ee[j]*x[i,j,s],j in JJ}) <= H[s]+Gt[i,s]*M[s]);
  @constraint(m,EndIndC[i in II, s in SS;s>=2],Gt[i,s] + Ft[i,s] == 1);
  @constraint(m,nonAT1[i in II, s in SS;s>=2],t[i,s] >= t[i,1] - (1 - F[i,s])*M[s]);
  @constraint(m,nonAT2[i in II, s in SS;s>=2],t[i,s] <= t[i,1] + (1 - F[i,s])*M[s]);
  @constraint(m,nonAX1[i in II, j in JJ, s in SS;s>=2],x[i,j,s] >= x[i,j,1] - (1 - F[i,s]));
  @constraint(m,nonAX2[i in II, j in JJ, s in SS;s>=2],x[i,j,s] <= x[i,j,1] + (1 - F[i,s]));
  @constraint(m,budget[s in SS],sum{x[i,j,s]*b[i,j] ,i in II, j in JJ} <= B);
  @constraint(m,oneAction[i in II, s in SS],sum{x[i,j,s],j in JJ} <= 1);
  @constraint(m,timecri1[i in II, s in SS;s>=2], F[i,s] + Gt[i,s] <= 1 + S1[i,s]);
  @constraint(m,timecri2[i in II, s in SS;s>=2], F[i,s] + Ft[i,s] <= 1 + S2[i,s]);
  @constraint(m,timecri3[i in II, s in SS;s>=2], G[i,s] + Gt[i,s] <= 1 + S3[i,s]);
  @constraint(m,timecriBind[i in II, s in SS;s>=2], S1[i,s] + S2[i,s] + S3[i,s] == 1);
  @constraint(m,linearZ11[i in II, j in JJ, s in SS;s>=2], Z1[i,j,s] <= S1[i,s]);
  @constraint(m,linearZ12[i in II, j in JJ, s in SS;s>=2], Z1[i,j,s] <= x[i,j,s]);
  @constraint(m,linearZ13[i in II, j in JJ, s in SS;s>=2], Z1[i,j,s] >= S1[i,s] + x[i,j,s] - 1);
  @constraint(m,linearZ21[i in II, j in JJ, s in SS;s>=2], Z2[i,j,s] <= S2[i,s]);
  @constraint(m,linearZ22[i in II, j in JJ, s in SS;s>=2], Z2[i,j,s] <= x[i,j,s]);
  @constraint(m,linearZ23[i in II, j in JJ, s in SS;s>=2], Z2[i,j,s] >= S2[i,s] + x[i,j,s] - 1);
  @constraint(m,linearZ31[i in II, j in JJ, s in SS;s>=2], Z3[i,j,s] <= S3[i,s]);
  @constraint(m,linearZ32[i in II, j in JJ, s in SS;s>=2], Z3[i,j,s] <= x[i,j,s]);
  @constraint(m,linearZ33[i in II, j in JJ, s in SS;s>=2], Z3[i,j,s] >= S3[i,s] + x[i,j,s] - 1);
  @constraint(m,linearY1[i in II, s in SS;s>=2], Y[i,s] <= M[s]*S1[i,s]);
  @constraint(m,linearY2[i in II, s in SS;s>=2], Y[i,s] <= t[i,s]);
  @constraint(m,linearY3[i in II, s in SS;s>=2], Y[i,s] >= t[i,s] - M[s]*(1-S1[i,s]));

  return m
end

function solveSub(xmaster,tmaster,D,rsub,Hsub,b,B,ee,II,JJ,GG,zint)
    # solve the lagrangian relaxation of the subproblem
    # categorize each activity and build the (relaxed) subproblem
    I1,I2,I3,tau = obtainIs(xmaster,tmaster,D,rsub,Hsub,b,B,ee,II,JJ,GG);

    # iteratively solve the lagrangian multiplier
    # initialization
    lambda1 = Dict{Any,Any}();
    lambda2 = Dict{Any,Any}();
    ppi = Dict{Any,Any}();
    for i in I1
        ppi[i] = 0;
        for j in JJ
            lambda1[i,j] = 0;
        end
    end
    for i in I2
        ppi[i] = 0;
        for j in JJ
            lambda2[i,j] = 0;
        end
    end
    # start the iteration
    stopBool = false;
    k = 0;
    μ = 2;
    zk = 0;
    counter = 0;
    mr = subLagP(xmaster,tmaster,D,Hsub,tau,b,B,ee,II,I1,I2,I3,JJ,GG,lambda1,lambda2,ppi);
    while (k<=20)&&(!stopBool)
        k += 1;
        @objective(mr,Min,mr.varDict[:tNr]+sum{lambda1[i,j]*(xmaster[i,j] - mr.varDict[:xr1][i,j]),i in I1,j in JJ}
            +sum{lambda2[i,j]*(xmaster[i,j]-mr.varDict[:xr2][i,j]),i in I2,j in JJ}
            +sum{ppi[i]*(tmaster[i]-mr.varDict[:ter][i]),i in I1}+sum{ppi[i]*(tmaster[i]-mr.varDict[:ter][i]),i in I2});
        solve(mr);
        xk1 = getvalue(mr.varDict[:xr1]);
        xk2 = getvalue(mr.varDict[:xr2]);
        tk = getvalue(mr.varDict[:ter]);
        if getobjectivevalue(mr) > zk
            counter = 0;
        else
            counter += 1;
        end
        if counter >= 5
            μ = μ/2;
        end
        zk = getobjectivevalue(mr);
        stopTemp = true;
        for j in JJ
            for i in I1
                if abs(xk1[i,j] - xmaster[i,j]) >= 1e-5
                    stopTemp = false;
                end
            end
            for i in I2
                if abs(xk2[i,j] - xmaster[i,j]) >= 1e-5
                    stopTemp = false;
                end
            end
        end
        for i in I1
            if abs(tk[i] - tmaster[i]) >= 1e-5
                stopTemp = false;
            end
        end
        for i in I2
            if abs(tk[i] - tmaster[i]) >= 1e-5
                stopTemp = false;
            end
        end
        if !stopTemp
            denomSum = 0;
            for i in I1
                denomSum += (tmaster[i] - tk[i])^2;
            end
            for i in I2
                denomSum += (tmaster[i] - tk[i])^2;
            end
            for j in JJ
                for i in I1
                    denomSum += (xmaster[i,j] - xk1[i,j])^2;
                end
                for i in I2
                    denomSum += (xmaster[i,j] - xk2[i,j])^2;
                end
            end
            # obtain the step length according to Held and Karp
            ν = (zint - zk)/denomSum*μ;
            # update the lambdas and ppis
            for i in I1
                ppi[i] += ν*(tmaster[i] - tk[i]);
            end
            for i in I2
                ppi[i] += ν*(tmaster[i] - tk[i]);
            end
            for j in JJ
                for i in I1
                    lambda1[i,j] += ν*(xmaster[i,j] - xk1[i,j]);
                end
                for i in I2
                    lambda2[i,j] += ν*(xmaster[i,j] - xk2[i,j]);
                end
            end
        else
            stopBool = true;
        end
    end
    lambda = merge(lambda1,lambda2);
    return ppi,lambda1,lambda2,lambda,zk;
end

function solveSub2(xmaster,tmaster,D,rsub,Hsub,b,B,ee,II,JJ,GG,zint,M)
    # solve the lagrangian relaxation of the subproblem
    # categorize each activity and build the (relaxed) subproblem
    I1,I2,I3,tau = obtainIs(xmaster,tmaster,D,rsub,Hsub,b,B,ee,II,JJ,GG);

    # iteratively solve the lagrangian multiplier
    # initialization
    lambda1 = Dict{Any,Any}();
    lambda2 = Dict{Any,Any}();
    ppi = Dict{Any,Any}();
    for i in I1
        ppi[i] = 0;
        for j in JJ
            lambda1[i,j] = 0;
        end
    end
    for i in I2
        ppi[i] = 0;
        for j in JJ
            lambda2[i,j] = 0;
        end
    end
    # start the iteration
    stopBool = false;
    k = 0;
    μ = 2;
    zk = 0;
    counter = 0;
    mr = subLag(xmaster,tmaster,D,rsub,Hsub,M[s],b,B,ee,II,I1,I2,I3,JJ,GG,lambda1,lambda2,ppi);
    while (k<=50)&&(!stopBool)
        k += 1;
        @objective(mr,Min,mr.varDict[:tNre]+sum{lambda1[i,j]*(xmaster[i,j] - mr.varDict[:xre1][i,j]),i in I1,j in JJ}
            +sum{lambda2[i,j]*(xmaster[i,j]-mr.varDict[:xre2][i,j]),i in I2,j in JJ}
            +sum{ppi[i]*(tmaster[i]-mr.varDict[:tere][i]),i in I1}+sum{ppi[i]*(tmaster[i]-mr.varDict[:tere][i]),i in I2});
        solve(mr);
        xk1 = getvalue(mr.varDict[:xre1]);
        xk2 = getvalue(mr.varDict[:xre2]);
        tk = getvalue(mr.varDict[:tere]);
        if getobjectivevalue(mr) > zk
            counter = 0;
        else
            counter += 1;
        end
        if counter >= 5
            μ = μ/2;
        end
        zk = getobjectivevalue(mr);
        stopTemp = true;
        for j in JJ
            for i in I1
                if abs(xk1[i,j] - xmaster[i,j]) >= 1e-5
                    stopTemp = false;
                end
            end
            for i in I2
                if abs(xk2[i,j] - xmaster[i,j]) >= 1e-5
                    stopTemp = false;
                end
            end
        end
        for i in I1
            if abs(tk[i] - tmaster[i]) >= 1e-5
                stopTemp = false;
            end
        end
        for i in I2
            if abs(tk[i] - tmaster[i]) >= 1e-5
                stopTemp = false;
            end
        end
        if !stopTemp
            denomSum = 0;
            for i in I1
                denomSum += (tmaster[i] - tk[i])^2;
            end
            for i in I2
                denomSum += (tmaster[i] - tk[i])^2;
            end
            for j in JJ
                for i in I1
                    denomSum += (xmaster[i,j] - xk1[i,j])^2;
                end
                for i in I2
                    denomSum += (xmaster[i,j] - xk2[i,j])^2;
                end
            end
            # obtain the step length according to Held and Karp
            ν = (zint - zk)/denomSum*μ;
            # update the lambdas and ppis
            for i in I1
                ppi[i] += ν*(tmaster[i] - tk[i]);
            end
            for i in I2
                ppi[i] += ν*(tmaster[i] - tk[i]);
            end
            for j in JJ
                for i in I1
                    lambda1[i,j] += ν*(xmaster[i,j] - xk1[i,j]);
                end
                for i in I2
                    lambda2[i,j] += ν*(xmaster[i,j] - xk2[i,j]);
                end
            end
        else
            stopBool = true;
        end
    end
    lambda = merge(lambda1,lambda2);
    return ppi,lambda1,lambda2,lambda,zk;
end
