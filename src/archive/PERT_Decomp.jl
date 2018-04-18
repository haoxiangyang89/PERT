# this is the code for decomposition of the stochastic program

using JuMP,Gurobi,Distributions;

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

function solveSub(xmaster,tmaster,D,rsub,Hsub,b,B,ee,II,JJ,GG)
    # solve the lagrangian relaxation of the subproblem
    # categorize each activity and build the (relaxed) subproblem
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

    # build the MIP subproblem
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
    solve(msub);
    zint = getobjectivevalue(msub);
    x1v = getvalue(x1);
    x2v = getvalue(x2);
    x3v = getvalue(x3);
    tv = getvalue(te);

    # build the Lagrangian Relaxation of the MIP subproblem
    mr = Model();
    @variables(mr,begin
        0 <= xr1[i in I1, j in JJ] <= 1
        0 <= xr2[i in I2, j in JJ] <= 1
        xr3[i in I3, j in JJ], Bin
        ter[i in II] >= 0
        tNr >= 0
    end);
#    @constraint(mr,binding1[i in I1],ter[i] == t[i]);
#    @constraint(mr,binding2[i in I2],ter[i] == t[i]);
    @constraint(mr,sepN[i in II], tNr >= ter[i]);
    @constraint(mr,sepi1[g in 1:length(GG);GG[g][1] in I1], ter[GG[g][2]] - ter[GG[g][1]] >= tau[GG[g][1]]);
    @constraint(mr,sepi2[g in 1:length(GG);GG[g][1] in I2], ter[GG[g][2]] - ter[GG[g][1]] >= tau[GG[g][1]]);
    @constraint(mr,sepi3[g in 1:length(GG);GG[g][1] in I3], ter[GG[g][2]] - ter[GG[g][1]] >= tau[GG[g][1]]*(1-sum{ee[j]*xr3[GG[g][1],j], j in JJ}));
    @constraint(mr,oneAction1[i in I1], sum{xr1[i,j], j in JJ} <= 1);
    @constraint(mr,oneAction2[i in I2], sum{xr2[i,j], j in JJ} <= 1);
    @constraint(mr,oneAction3[i in I3], sum{xr3[i,j], j in JJ} <= 1);
    @constraint(mr,budget,sum{xr1[i,j]*b[i,j], i in I1, j in JJ} +
      sum{xr2[i,j]*b[i,j], i in I2, j in JJ} + sum{xr3[i,j]*b[i,j], i in I3, j in JJ}<= B);

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
    while (k<=50)&&(!stopBool)
        k += 1;
        @objective(mr,Min,tNr+sum{lambda1[i,j]*(xmaster[i,j] - xr1[i,j]),i in I1,j in JJ}
            +sum{lambda2[i,j]*(xmaster[i,j]-xr2[i,j]),i in I2,j in JJ}
            +sum{ppi[i]*(tmaster[i]-ter[i]),i in I1}+sum{ppi[i]*(tmaster[i]-ter[i]),i in I2});
#        @expression(mr,lagrangianD,tNr+sum{lambda1[i,j]*(-xr1[i,j]),i in I1,j in JJ}
#            +sum{lambda2[i,j]*(-xr2[i,j]),i in I2,j in JJ});
        solve(mr);
        xk1 = getvalue(xr1);
        xk2 = getvalue(xr2);
        tk = getvalue(ter);
#        zD = getvalue(lagrangianD);
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
            ν = (zint*1.01 - zk)/denomSum*μ;
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

    # obtain the optimal lambda1 and lambda2
    # using those lambdas to generate a linear cut for t
#    mrr = Model();
#    @variables(mrr,begin
#        0 <= xrr1[i in I1, j in JJ] <= 1
#        0 <= xrr2[i in I2, j in JJ] <= 1
#        0 <= xrr3[i in I3, j in JJ] <= 1
#        terr[i in II] >= 0
#        tNrr >= 0
#    end);
#    @constraint(mrr,bindingrr1[i in I1],terr[i] == t[i]);
#    @constraint(mrr,bindingrr2[i in I2],terr[i] == t[i]);
#    @constraint(mrr,sepN[i in II], tNrr >= terr[i]);
#    @constraint(mrr,sepi1[g in 1:length(GG);GG[g][1] in I1], terr[GG[g][2]] - terr[GG[g][1]] >= tau[GG[g][1]]);
#    @constraint(mrr,sepi2[g in 1:length(GG);GG[g][1] in I2], terr[GG[g][2]] - terr[GG[g][1]] >= tau[GG[g][1]]);
#    @constraint(mrr,sepi3[g in 1:length(GG);GG[g][1] in I3], terr[GG[g][2]] - terr[GG[g][1]] >= tau[GG[g][1]]*(1-sum{ee[j]*xrr3[GG[g][1],j], j in JJ}));
#    @constraint(mrr,oneAction1[i in I1], sum{xrr1[i,j], j in JJ} <= 1);
#    @constraint(mrr,oneAction2[i in I2], sum{xrr2[i,j], j in JJ} <= 1);
#    @constraint(mrr,oneAction3[i in I3], sum{xrr3[i,j], j in JJ} <= 1);
#    @constraint(mrr,budget,sum{xrr1[i,j]*b[i,j], i in I1, j in JJ} +
#      sum{xrr2[i,j]*b[i,j], i in I2, j in JJ} + sum{xrr3[i,j]*b[i,j], i in I3, j in JJ}<= B);
#    @objective(mrr,Min,tNrr+sum{lambda1[i,j]*(x[i,j] - xrr1[i,j]),i in I1,j in JJ}
#        +sum{lambda2[i,j]*(x[i,j]-xrr2[i,j]),i in I2,j in JJ});
#    solve(mrr);
#    zlin = getobjectivevalue(mrr);
#    pi = Dict{Any,Any}();
#    for i in I1
#        pi[i] = getdual(bindingrr1[i]);
#    end
#    for i in I2
#        pi[i] = getdual(bindingrr2[i]);
#    end
    lambda = merge(lambda1,lambda2);
    return ppi,lambda,zint,zk,x1v,x2v,x3v,tv;
end

function solveMaster(InputAdd)
    D,r,H,b,B,ee,II,JJ,M,SS,GG,dH,dR,p = readIn(InputAdd);
    # set up the solver as Gurobi
    env = Gurobi.Env();
    setparams!(env;MIPGap = 1e-8,MIPFocus = 2);
    solver = GurobiSolver(env);

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

    # Step 0: initialization
    k = 1;
    l = zeros(length(SS));
    ub = Any[9999999];
    lb = Any[0];
    keepIter = true;
    ϵ = 0.01;
    counter = 0;
    maxIter = 5;
    cuts = Dict{Any,Any}();
    xrec = [];

    while keepIter
        # Step 1: solve the master and obtain the incumbent solutions
        solve(m);
        xhat = getvalue(x)[:,:];
        push!(xrec,xhat);
        that = getvalue(t)[:];
        fv = getobjectivevalue(m);
        recourseFV = Dict{Any,Any}();
        recourseFV[1] = getvalue(tN);
        recInt = Dict{Any,Any}();
        recInt[1] = getvalue(tN);
        # update the lower bound
        if fv > lb[length(lb)]
            push!(lb,fv);
        else
            temp = lb[length(lb)];
            push!(lb,temp);
        end
        # Step 2: for each scenario
        # solve the subproblem
        for s in SS[2:length(SS)]
            # generate the cut
            rscen = Dict{Any,Any}();
            for kr in keys(r)
                if kr[2] == s
                    rscen[kr[1]] = r[kr];
                end
            end
            ppi,lambda,intFV,rFV,x1s,x2s,x3s,ts = solveSub(xhat,that,D,rscen,H[s],b,B,ee,II,JJ,GG);
            recInt[s] = intFV;
            recourseFV[s] = rFV;
            # append the cut to the master
            l[s] += 1;
            cuts[s,l[s]] = @constraint(m,θ[s] >= sum{ppi[i]*(t[i] - that[i]),i in II;i in keys(ppi)}
                + sum{lambda[i,j]*(x[i,j] - xhat[i,j]), i in II, j in JJ;(i,j) in keys(lambda)} + rFV);
            #@constraint(m,θ[s] >= sum{lambda[i,j]*(x[i,j] - xhat[i,j]), i in II, j in JJ;(i,j) in keys(lambda)} + rFV);
        end
        # Step 3: update the upper bound
        ubtemp = p*recInt[1] + (1-p)/(length(SS)-1)*sum([recInt[s] for s in 2:length(SS)]);
        if ubtemp < ub[length(ub)]
          push!(ub,ubtemp);
          xfin = xhat;
          tfin = that;
        else
          push!(ub,ub[length(ub)]);
        end
        # check if it meets the criterion to stop
        k += 1;
        if ub[k] - lb[k] >= ϵ
            if (ub[k] - lb[k]) <= (ub[k-1] - lb[k-1]) - ϵ
                counter = 0;
            else
                counter += 1;
            end
            if counter >= maxIter
                keepIter = false;
            end
        else
            keepIter = false;
        end
    end
    # record the final solutions
    fvfin = ub[length(ub)];
    return xfin,tfin,fvfin;
end

#xf,tf,fvf = solveMaster("test_Input_graph_Full.csv");
#print(xf,tf,fvf);
