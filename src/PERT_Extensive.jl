# This is the code for stochastic program with extensive formulation

using JuMP,Gurobi,Distributions;

# readin the data
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
    elseif Hname == "Normal"
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

# generate optimal solution of full-stochastc case
function fullStoch(InputAdd)
    D,r,H,b,B,ee,II,JJ,M,SS,GG,dH,dR,p = readIn(InputAdd);
    m = Model();
    # set up the solver as Gurobi
    solver = GurobiSolver();

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
    @constraint(m,sepi[g in 1:length(GG),s in SS],
      t[GG[g][2],s] - t[GG[g][1],s] >=
      (r[GG[g][1],s] - 1)*Y[GG[g][1],s] + (r[GG[g][1],s]*D[GG[g][1]]+H[s]*(1-r[GG[g][1],s]))*S1[GG[g][1],s]
      - (r[GG[g][1],s]*D[GG[g][1]]) * sum{(ee[j]*Z1[GG[g][1],j,s]),j in JJ}
      + D[GG[g][1]]*(S2[GG[g][1],s] - sum{ee[j]*Z2[GG[g][1],j,s],j in JJ})
      + D[GG[g][1]]*r[GG[g][1],s]*(S3[GG[g][1],s] - sum{ee[j]*Z3[GG[g][1],j,s],j in JJ})
      );
    @constraint(m,sepN[i in II, s in SS],tN[s] >= t[i,s]);
    @constraint(m,timeInd1[i in II, s in SS],t[i,s] >= H[s]-F[i,s]*H[s]);
    @constraint(m,timeInd2[i in II, s in SS],t[i,s] <= H[s]+G[i,s]*M[s]);
    @constraint(m,timeIndC[i in II, s in SS],G[i,s] + F[i,s] == 1);
    @constraint(m,EndInd1[i in II, s in SS],t[i,s] + D[i]*(1-sum{ee[j]*x[i,j,s],j in JJ}) >= H[s]-Ft[i,s]*H[s]);
    @constraint(m,EndInd2[i in II, s in SS],t[i,s] + D[i]*(1-sum{ee[j]*x[i,j,s],j in JJ}) <= H[s]+Gt[i,s]*M[s]);
    @constraint(m,EndIndC[i in II, s in SS],Gt[i,s] + Ft[i,s] == 1);
    @constraint(m,nonAT1[i in II, s in SS],t[i,s] >= t[i,1] - (1 - F[i,s])*M[s]);
    @constraint(m,nonAT2[i in II, s in SS],t[i,s] <= t[i,1] + (1 - F[i,s])*M[s]);
    @constraint(m,nonAX1[i in II, j in JJ, s in SS],x[i,j,s] >= x[i,j,1] - (1 - F[i,s]));
    @constraint(m,nonAX2[i in II, j in JJ, s in SS],x[i,j,s] <= x[i,j,1] + (1 - F[i,s]));
    @constraint(m,budget[s in SS],sum{x[i,j,s]*b[i,j] ,i in II, j in JJ} <= B);
    @constraint(m,oneAction[i in II, s in SS],sum{x[i,j,s],j in JJ} <= 1);
    @constraint(m,timecri1[i in II, s in SS], F[i,s] + Gt[i,s] <= 1 + S1[i,s]);
    @constraint(m,timecri2[i in II, s in SS], F[i,s] + Ft[i,s] <= 1 + S2[i,s]);
    @constraint(m,timecri3[i in II, s in SS], G[i,s] + Gt[i,s] <= 1 + S3[i,s]);
    @constraint(m,timecriBind[i in II, s in SS], S1[i,s] + S2[i,s] + S3[i,s] == 1);
    @constraint(m,linearZ11[i in II, j in JJ, s in SS], Z1[i,j,s] <= S1[i,s]);
    @constraint(m,linearZ12[i in II, j in JJ, s in SS], Z1[i,j,s] <= x[i,j,s]);
    @constraint(m,linearZ13[i in II, j in JJ, s in SS], Z1[i,j,s] >= S1[i,s] + x[i,j,s] - 1);
    @constraint(m,linearZ21[i in II, j in JJ, s in SS], Z2[i,j,s] <= S2[i,s]);
    @constraint(m,linearZ22[i in II, j in JJ, s in SS], Z2[i,j,s] <= x[i,j,s]);
    @constraint(m,linearZ23[i in II, j in JJ, s in SS], Z2[i,j,s] >= S2[i,s] + x[i,j,s] - 1);
    @constraint(m,linearZ31[i in II, j in JJ, s in SS], Z3[i,j,s] <= S3[i,s]);
    @constraint(m,linearZ32[i in II, j in JJ, s in SS], Z3[i,j,s] <= x[i,j,s]);
    @constraint(m,linearZ33[i in II, j in JJ, s in SS], Z3[i,j,s] >= S3[i,s] + x[i,j,s] - 1);
    @constraint(m,linearY1[i in II, s in SS], Y[i,s] <= M[s]*S1[i,s]);
    @constraint(m,linearY2[i in II, s in SS], Y[i,s] <= t[i,s]);
    @constraint(m,linearY3[i in II, s in SS], Y[i,s] >= t[i,s] - M[s]*(1-S1[i,s]));

    solve(m);
    objV = getobjectivevalue(m);
    xVal = getvalue(x);
    tVal = getvalue(t);
    return objV,tVal[:,1],xVal[:,:,1],H;
end

# generate optimal solution of deterministic case without disruption
function deterOpt(InputAdd)
    D,r,H,b,B,ee,II,JJ,M,SS,GG,dH,dR,p = readIn(InputAdd);
    m = Model();
    # set up the solver as Gurobi
    solver = GurobiSolver();

    @variables(m,begin
        t[i in II] >= 0
        tN >= 0
        x[i in II, j in JJ], Bin
    end);
    @objective(m,Min,tN);
    @constraint(m,sepN[i in II],tN>=t[i]);
    @constraint(m,sepi[g in 1:length(GG)],t[GG[g][2]]-t[GG[g][1]]>=D[GG[g][1]]*(1-sum{ee[j]*x[GG[g][1],j],j in JJ}));
    @constraint(m,oneAction[i in II],sum{x[i,j],j in JJ} <= 1);
    @constraint(m,budget,sum{x[i,j]*b[i,j],i in II, j in JJ} <= B);

    solve(m);
    optV = getobjectivevalue(m);
    tV = getvalue(t);
    xV = getvalue(x);
    return optV,tV[:],xV[:,:];
end

# evaluate the solution x coupled with nominal first stage t
function evalSub(x,t,N,B,ee,D,HList,rEval,II,JJ,GG)
    fVal = zeros(N);
    solver = GurobiSolver();
    for n in 1:N
        m = Model();
        I1 = [];
        I2 = [];
        I3 = [];
        tau = Dict{Any,Any}();
        for i in II
            if t[i] >= HList[n]
                push!(I3,i);
                tau[i] = D[i]*rEval[i][n];
            else
                if t[i] + (1-sum([ee[j]*x[i,j] for j in JJ]))*D[i] <= HList[n]
                    tau[i] = D[i]*(1-sum([ee[j]*x[i,j] for j in JJ]));
                    push!(I2,i);
                else
                    push!(I1,i);
                    tau[i] = (HList[n]-t[i]) + rEval[i][n]*(t[i]+D[i]*(1-sum([ee[j]*x[i,j] for j in JJ]))-HList[n]);
                end
            end
        end
        @variables(m,begin
            te[i in II] >= 0
            tN >= 0
            xe[i in II,j in JJ], Bin
        end);
        @objective(m,Min,tN);
        @constraint(m,bindingT1[i in I1],te[i] >= t[i]);
        @constraint(m,bindingT2[i in I2],te[i] >= t[i]);
        @constraint(m,bindingX1[i in I1,j in JJ],xe[i,j] >= x[i,j]);
        @constraint(m,bindingX2[i in I2,j in JJ],xe[i,j] >= x[i,j]);
        @constraint(m,sep1[g in GG;g[1] in I1],te[g[2]] - te[g[1]] >= tau[g[1]]);
        @constraint(m,sep2[g in GG;g[1] in I2],te[g[2]] - te[g[1]] >= tau[g[1]]);
        @constraint(m,sep3[g in GG;g[1] in I3],te[g[2]] - te[g[1]] >= tau[g[1]]*(1-sum([ee[j]*xe[g[1],j] for j in JJ])));
        @constraint(m,oneAction[i in II],sum{xe[i,j], j in JJ} <= 1);
        @constraint(m,budget,sum{xe[i,j],i in II, j in JJ} <= B);
        @constraint(m,sepN[i in II], tN >= te[i]);
        solve(m);
        fVal[n] = getobjectivevalue(m);
    end
    return fVal;
end

# generate the scenario for evaluation
function genScen(N,InputAdd)
    D,r,H,b,B,ee,II,JJ,M,SS,GG,dH,dR,p = readIn(InputAdd);
    # generate a list of H
    HList = rand(dH,N);
    # generate the magnitude of the disruption
    rEval = Dict{Any,Any}(0 => zeros(N));
    for i in II
        rEval[i] = round(rand(dR[i],N),3);
    end
    return N,B,ee,D,HList,rEval,II,JJ,GG,p;
end

function runExpr(N,InputAddF,InputAddS1,InputAddS2,InputAddD,option)
    N,B,ee,D,HList,rEval,II,JJ,GG,p = genScen(N,InputAddF);
    HList0 = [99999];
    rEval0 = Dict{Any,Any}(0 => zeros(1));
    for i in II
        rEval0[i] = [1];
    end

    # fully stochastic: both timing and magnitude stochastic
    @time ov1,tv1,xv1,H = fullStoch(InputAddF);
    f1 = evalSub(xv1,tv1,N,B,ee,D,HList,rEval,II,JJ,GG);
    fnoD1 = evalSub(xv1,tv1,1,B,ee,D,HList0,rEval0,II,JJ,GG);
    # semi-stochastic 1: only magnitude stochastic
    @time ov2,tv2,xv2,H = fullStoch(InputAddS1);
    f2 = evalSub(xv2,tv2,N,B,ee,D,HList,rEval,II,JJ,GG);
    fnoD2 = evalSub(xv2,tv2,1,B,ee,D,HList0,rEval0,II,JJ,GG);
    # semi-stochastic 2: only timing stochastic
    @time ov3,tv3,xv3,H = fullStoch(InputAddS2);
    f3 = evalSub(xv3,tv3,N,B,ee,D,HList,rEval,II,JJ,GG);
    fnoD3 = evalSub(xv3,tv3,1,B,ee,D,HList0,rEval0,II,JJ,GG);
    # semi-stochastic 1: both timing/magnitude deterministic, but with disruption
    @time ov4,tv4,xv4,H = fullStoch(InputAddD);
    f4 = evalSub(xv4,tv4,N,B,ee,D,HList,rEval,II,JJ,GG);
    fnoD4 = evalSub(xv4,tv4,1,B,ee,D,HList0,rEval0,II,JJ,GG);
    # deterministic case: no disruption
    @time ov5,tv5,xv5 = deterOpt(InputAddF);
    f5 = evalSub(xv5,tv5,N,B,ee,D,HList,rEval,II,JJ,GG);
    fnoD5 = evalSub(xv5,tv5,1,B,ee,D,HList0,rEval0,II,JJ,GG);
    return f1,fnoD1,f2,fnoD2,f3,fnoD3,f4,fnoD4,f5,fnoD5;
end
