# This is a collections of all functions
using JuMP,Distributions,CPLEX;

function readIn(InputAdd)
    # no of activities
    data = readdlm(InputAdd,',');
    # number of events
    lI = data[1,1];
    II = 1:lI;
    # number of crashing options
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
    d = Dict();
    H = Dict();
    # no disruption in first scenario

    # first read in the distribution names and then the parameters
    # information about the disruption time
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

    # information about the disruption magnitude
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
            d[i,1] = 0;
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
            d[i,1] = 0;
            disRList[i] = LogNormal(μ,σ);
        end
    end

    # the number of scenarios
    row += 1;
    S = data[row,1] + 1;
    # the probability of no disruption
    p = data[row,2];
    M = Dict();

    for s in 2:S
      # generate a disruption time
        H[s] = round(rand(distrH),4);
        # if the disruption time is larger than the entire project span, resample
        while H[s] > 97
            H[s] = round(rand(distrH),4);
        end
        for i in II
            d[i,s] = round(rand(disRList[i]),4);
        end
        M[s] = sum(values(D))*maximum([d[i,s] for i in II]);
    end
    H[1] = maximum(values(M));
    M[1] = 2*sum(values(D));
    # set of scenarios
    SS = 1:S;

    return D,d,H,b,B,ee,II,J,M,SS,G,distrH,disRList,p
end

# function to create the basic version of the master program
function createMaster(D,d,b,B,ee,II,JJ,SS,GG,p)
  mp = Model(solver = CplexSolver(CPXPARAM_ScreenOutput = 0));
  @variables(mp, begin
    θ[s in SS;s >= 2]
    x[i in II,j in JJ], Bin
    t[i in II] >= 0
    tN >= 0
  end);
  @constraint(mp, budgetConstr, sum(sum(b[i,j]*x[i,j] for j in JJ) for i in II) <= B);
  @constraint(mp, durationConstr[g in GG], t[g[2]] - t[g[1]] >= D[g[1]]*(1-sum(ee[j]*x[g[1],j] for j in JJ)));
  @constraint(mp, tNConstr[i in II], tN >= t[i]);
  @constraint(mp, xConstr[i in II], sum(x[i,j] for j in JJ) <= 1);
  @objective(mp,Min,p*tN + sum(1/(length(SS)-1)*θ[s] for s in SS if s >= 2));

  return mp;
end

# function to append Lagrangian cuts to the given program
function appendLCcuts(mp,cI,II,JJ)
  # cI(constraint item) contains πl, λl and constant
  @constraint(mp,mp.varDict[:θ][cI.s] >= sum(cI.πl[i]*mp.varDict[:t][i] for i in II if i in keys(πl))
                  + sum(sum(cI.λl[i,j]*mp.varDict[:x][i,j] for j in JJ if (i,j) in keys(λl)) for i in II)
                  + cI.z);
  return mp;
end

# function to append branch-and-bound cuts to the given program
function appendBNBcuts(mp,i,td,sign)
  if sign == 1
    @constraint(mp,mp.varDict[:t][i] <= td - 1e-4);
  else
    @constraint(mp,mp.varDict[:t][i] >= td);
  end
  return mp;
end

# function to obtain the status of each activity: before/after the disruption
function obtainIs(xs,ts,td,II)
  I1 = [];
  I2 = [];
  # if ts[i] < td, then i in I1, else i in I2
  for i in II
    if ts[i] < td
      push!(I1,i);
    else
      push!(I2,i);
    end
  end
  return I1,I2;
end

# function to build the MIP sub problem
function subInt(xs,ts,D,b,B,ee,II,I1,I2,dscen,JJ,GG)
  sip = Model(solver = CplexSolver());
  @variables(sip, begin
    x[i in II,j in JJ], Bin
    t[i in II] >= 0
    tN >= 0
  end);
  @constraint(mp, budgetConstr, sum(sum(b[i,j]*x[i,j] for j in JJ) for i in II) <= B);
  @constraint(mp, durationConstr1[g in GG;g[1] in I1], t[g[2]] - t[g[1]] >= D[g[1]]*(1-sum(ee[j]*x[g[1],j] for j in JJ)));
  @constraint(mp, durationConstr2[g in GG;g[1] in I2], t[g[2]] - t[g[1]] >= (D[g[1]]+dscen[g[1]])*(1-sum(ee[j]*x[g[1],j] for j in JJ)));
  @constraint(mp, tNConstr[i in II], tN >= t[i]);
  @constraint(mp, xConstr[i in II], sum(x[i,j] for j in JJ) <= 1);
  @constraint(mp, tConstant[i in I1], t[i] == ts[i]);
  @constraint(mp, xConstant[i in I1, j in JJ], x[i,j] == xs[i,j]);
  return sip;
end

# function to solve for the coefficents used in Lagrangian cuts
# the easier version of the sub problem: without big M
function subLagP(D,dscen,b,B,ee,II,I1,I2,JJ,GG)
  spLp = Model(solver = CplexSolver());
  @variables(spLp, begin
    tN >= 0
    t[i in II] >= 0
    0 <= x1[i in I1, j in JJ] <= 1
    x2[i in I2, j in JJ], Bin
  end);
  @constraint(spLp,tNConstr[i in II],tN >= t[i]);
  @constraint(spLp,durationConstr1[g in GG; g[1] in I1], t[g[2]] - t[g[1]] >= D[g[1]]*(1-sum(ee[j]*x1[g[1],j] for j in JJ)));
  @constraint(spLp,durationConstr2[g in GG; g[1] in I2], t[g[2]] - t[g[1]] >= (D[g[1]]+dscen[g[1]])*(1-sum(ee[j]*x2[g[1],j] for j in JJ)));
  @constraint(spLp,xConstr1[i in I1], sum(x1[i,j] for j in JJ) <= 1);
  @constraint(spLp,xConstr2[i in I1], sum(x2[i,j] for j in JJ) <= 1);
  @constraint(spLp,budgetConstr, sum(sum(b[i,j]*x1[i,j] for j in JJ) for i in I1) + sum(sum(b[i,j]*x2[i,j] for j in JJ) for i in I2) <= B);

  return spLp;
end

# function to solve the Lagrangian relaxation with big M
function subLag(xs,ts,D,dscen,td,M,b,B,ee,II,I1,I2,JJ,GG,πle,λle)
  spL = Model(solver = CplexSolver());
  @variables(spL, begin
    tN >= 0
    # t and x are the auxiliary variables with the added equality constraints
    # τ and z are the actual subproblem variables
    t[i in II] >= 0
    τ[i in II] >= 0
    z[i in II], Bin
    0 <= x[i in II] <= 1
    F[i in II], Bin
    G[i in II], Bin
    # S is the linearization of the G[i]*z[i,j]
    S[i in II, j in JJ], Bin
  );
  @constraint(spL, tNConstr[i in II], tN >= τ[i]);
  @constraint(spL, FLogic[i in II], td - F[i]*M <= t[i]);
  @constraint(spL, GLogic[i in II], td + G[i]*M >= t[i]);
  @constraint(spL, FGLogic[i in II], F[i] + G[i] == 1);
  @constraint(spL, τLogic1[i in II], τ[i] + (1 - F[i])*M >= t[i]);
  @constraint(spL, τLogic2[i in II], τ[i] - (1 - F[i])*M <= t[i]);
  @constraint(spL, zlogic1[i in II, j in JJ], z[i,j] + 1 - F[i] >= x[i,j]);
  @constraint(spL, zlogic2[i in II, j in JJ], z[i,j] - 1 + F[i] <= x[i,j]);
  @constraint(spL, durationConstr[g in GG], τ[g[2]] - τ[g[1]] >= (D[g[1]] + dscen[g[1]]*G[g[1]]) -
                    sum(D[g[1]]*ee[j]*z[g[1],j] for j in JJ) - sum(dscen[g[1]]*ee[j]*S[g[1],j] for j in JJ));
  @constraint(spL, xConstr[j in JJ], sum(z[i,j] for j in JJ) <= 1);
  @constraint(spL, budgetConstr, sum(sum(b[i,j]*z[i,j] for j in JJ) for i in II));
  # linearization constraints
  @constraint(spL, Slinear1[i in II, j in JJ], S[i,j] <= G[i]);
  @constraint(spL, Slinear2[i in II, j in JJ], S[i,j] <= z[i,j]);
  @constraint(spL, Slinear3[i in II, j in JJ], S[i,j] >= G[i] + z[i,j] - 1);

  @objective(spL, Min, tN - sum(πle[i]*(t[i] - ts[i]) for i in I1) - sum(sum(λle[i,j]*(x[i,j] - xs[i,j]) for j in JJ) for i in I1));
  return spL;
end

# solve the lagrangian relaxation of the subproblem
# categorize each activity and build the (relaxed) subproblem
function solveSub(xs,ts,D,dscen,td,b,B,ee,II,JJ,GG,zint)
  I1,I2 = obtainIs(xs,ts,td,II);

  # iteratively solve the lagrangian multiplier
  # initialization with every Lagrangian multiplier as 0
  λls = Dict();
  πls = Dict();
  for i in I1
    πls[i] = 0;
    for j in JJ
      λls[i,j] = 0;
    end
  end

  # start the iteration
  stopBool = false;
  k = 0;
  μ = 2;
  zk = 0;
  counter = 0;
  mr = subLagP(xs,ts,D,dscen,td,b,B,ee,II,I1,I2,JJ,GG);
  while (k<=50)&(!stopBool)
    k += 1;
    @objective(mr,Min,mr.varDict[:tN] - sum(πls[i]*(mr.varDict[:t][i] - ts[i]) for i in I1)
                - sum(sum(λls[i,j]*(mr.varDict[:x1][i,j] - xs[i,j]) for j in JJ) for i in I1));
    solve(mr);
    xk1 = getvalue(mr.varDict[:x1]);
    tk = getvalue(mr.varDict[:t]);
    if getobjectivevalue(mr) > zk
        counter = 0;
    else
        counter += 1;
    end
    if counter >= 5
        μ = μ/2;
    end
    zk = getobjectivevalue(mr);

    # check if there is still room for improvement
    stopTemp = true;
    for j in JJ
      for i in I1
        if abs(xk1[i,j] - xs[i,j]) >= 1e-5
            stopTemp = false;
        end
      end
    end
    for i in I1
      if abs(tk[i] - ts[i]) >= 1e-5
        stopTemp = false;
      end
    end
    # if there is still room for improvement, calculate the direction of improving λ and π
    if !stopTemp
        denomSum = 0;
        for i in I1
          denomSum += (ts[i] - tk[i])^2;
        end
        for j in JJ
          for i in I1
            denomSum += (xs[i,j] - xk1[i,j])^2;
          end
        end
        # obtain the step length according to Held and Karp
        ν = (zint - zk)/denomSum*μ;
        # update the lambdas and ppis
        for i in I1
            πls[i] += ν*(ts[i] - tk[i]);
        end
        for j in JJ
          for i in I1
            λls[i,j] += ν*(xs[i,j] - xk1[i,j]);
          end
        end
    else
      stopBool = true;
    end
  end

  # now we have a set of Lagrangian coefficients
  @objective(mr,Min,mr.varDict[:tN] - sum(πls[i]*(mr.varDict[:t][i] - ts[i]) for i in I1)
              - sum(sum(λls[i,j]*(mr.varDict[:x1][i,j] - xs[i,j]) for j in JJ) for i in I1));
  solve(mr);
  zlags = getobjectivevalue(mr);
  return πls,λls,zlags;
end

# function to generate a Lagrangian cut
function generateCut(s,πg,λg,zlag,xSol,tSol,II,JJ)
  lc = LagCut(s,πg,λg,zlag - sum(πg[i]*tSol[i] for i in II) - sum(sum(λg[i,j]*xSol[i,j] for j in JJ) for i in II));
  return lc;
end
