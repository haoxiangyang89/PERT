# dumpster
# wrong evaluation: nonbinding before the disruption
function evalXH_1(x,N,InputAdd)
    D,r,H,b,B,ee,II,JJ,M,SS,GG,μ,p,la,ub = readIn(InputAdd);
    fVal = zeros(N);
    distrH = Exponential(μ);
    # generate a list of H
    HList = rand(distrH,N);
    Itemp = Array(II);
    for g in GG
        if g[2] in Itemp
            deleteat!(Itemp,findin(Itemp,[g[2]]));
        end
    end
    # generate the magnitude of the disruption
    rEval = Dict{Any,Any}(0 => zeros(N));
    for i in II
        distrAB = Uniform(la[i],ub[i]+1e-10);
        rEval[i] = round(rand(distrAB,N),3);
    end
    for n in 1:N
        # add a source node
        for item in Itemp
            push!(GG,[0 item]);
        end
        # initialize with a max time
        t = Dict{Any,Any}(0 => 0);
        D[0] = 0;
        for i in II
            t[i] = Float64(1.0);
        end
        # start the Dijkstra
        cN = 0;
        histList = [];
        # we assume in the input the last node is the sink
        while length(histList)<length(II)+1
            # add the current node to the list which has been searched
            push!(histList,cN);
            # search the current node
            if cN == 0
                tau = 0;
            else
                if -t[cN] >= HList[n]
                    tau = D[cN]*(1-sum([ee[j]*x[cN,j] for j in JJ]))*rEval[cN][n];
                elseif -t[cN] + (1-sum([ee[j]*x[cN,j] for j in JJ]))*D[cN] <= HList[n]
                    tau = D[cN]*(1-sum([ee[j]*x[cN,j] for j in JJ]));
                else
                    tau = (HList[n]+t[cN]) + rEval[cN][n]*(-t[cN]+D[cN]*(1-sum([ee[j]*x[cN,j] for j in JJ]))-HList[n]);
                end
            end
            for g in GG
                if g[1] == cN
                    if t[g[2]] > t[g[1]] - tau
                        t[g[2]] = t[g[1]] - tau;
                    end
                end
            end
            # search for the smallest t
            minI = 1;
            minInd = 99;
            for i in II
                if !(i in histList)
                    if minI > t[i]
                        minI = t[i];
                        minInd = i;
                    end
                end
            end
            cN = minInd;
        end
        fVal[n] = -t[length(II)];
    end
    return fVal;
end

# generate optimal solution of semi-stochastic case
function semiStoch1(InputAdd)
    D,r,H,b,B,ee,II,JJ,M,SS,GG,μ,p,la,ub = readIn(InputAdd);
    m = Model();
    # set up the solver as Gurobi
    solver = GurobiSolver();

    # set up the variables
    SS = [1,2];
    H[2] = μ*(1-p);
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

    @objective(m,Min,sum{tN[s], s in SS;s>1});
    @constraint(m,sepi[g in 1:length(GG),s in SS],
      t[GG[g][2],s] - t[GG[g][1],s] >= (r[GG[g][1],s] - 1)*Y[GG[g][1],s] + (r[GG[g][1],s]*D[GG[g][1]]+H[s]*(1-r[GG[g][1],s]))*S1[GG[g][1],s]
      - sum{(r[GG[g][1],s]*D[GG[g][1]])*ee[j]*Z1[GG[g][1],j,s],j in JJ}
      + D[GG[g][1]]*(S2[GG[g][1],s] - sum{ee[j]*Z2[GG[g][1],j,s],j in JJ}) + D[GG[g][1]]*r[GG[g][1],s]*(S3[GG[g][1],s] - sum{ee[j]*Z3[GG[g][1],j,s],j in JJ})
      );
    @constraint(m,sepN[i in II, s in SS],tN[s] >= t[i,s]);
    @constraint(m,timeInd1[i in II, s in SS],t[i,s] >= H[s]-F[i,s]*M);
    @constraint(m,timeInd2[i in II, s in SS],t[i,s] <= H[s]+G[i,s]*M);
    @constraint(m,timeIndC[i in II, s in SS],G[i,s] + F[i,s] == 1);
    @constraint(m,EndInd1[i in II, s in SS],t[i,s] + D[i]*(1-sum{ee[j]*x[i,j,s],j in JJ}) >= H[s]-Ft[i,s]*M);
    @constraint(m,EndInd2[i in II, s in SS],t[i,s] + D[i]*(1-sum{ee[j]*x[i,j,s],j in JJ}) <= H[s]+Gt[i,s]*M);
    @constraint(m,EndIndC[i in II, s in SS],Gt[i,s] + Ft[i,s] == 1);
    @constraint(m,budget[s in SS],sum{x[i,j,s]*b[i,j], i in II, j in JJ} <= B);
    @constraint(m,oneAction[i in II, s in SS],sum{x[i,j,s],j in JJ} <= 1);
    @constraint(m,nonAT1[i in II, s in SS;s>1],t[i,s] >= t[i,2] - (1 - F[i,s])*M);
    @constraint(m,nonAT2[i in II, s in SS;s>1],t[i,s] <= t[i,2] + (1 - F[i,s])*M);
    @constraint(m,nonAX1[i in II, j in JJ, s in SS;s>1],x[i,j,s] >= x[i,j,2] - (1 - F[i,s]));
    @constraint(m,nonAX2[i in II, j in JJ, s in SS;s>1],x[i,j,s] <= x[i,j,2] + (1 - F[i,s]));
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
    @constraint(m,linearY1[i in II, s in SS], Y[i,s] <= M*S1[i,s]);
    @constraint(m,linearY2[i in II, s in SS], Y[i,s] <= t[i,s]);
    @constraint(m,linearY3[i in II, s in SS], Y[i,s] >= t[i,s] - M*(1-S1[i,s]));

    solve(m);
    objV = getobjectivevalue(m);
    xVal = getvalue(x);
    tVal = getvalue(t);
    return objV,tVal[:,2],xVal[:,:,2],H;
end

# given x, generate the non-disruption scenario t
function evalT(x,InputAdd)
    D,r,H,b,B,ee,II,JJ,M,SS,GG,μ,p,la,ub = readIn(InputAdd);
    Itemp = Array(II);
    for g in GG
        if g[2] in Itemp
            deleteat!(Itemp,findin(Itemp,[g[2]]));
        end
    end

    for item in Itemp
        push!(GG,[0 item]);
    end
    # initialize with a max time
    t = Dict{Any,Any}(0 => 0);
    D[0] = 0;
    for i in II
        t[i] = Float64(1.0);
    end
    # start the Bellman-Ford
    # we assume in the input the last node is the sink
    # also assume no disruption here
    for i in length(II)
        for g in GG
            if g[1] != 0
                if t[g[1]] - D[g[1]]*(1 - sum([ee[j]*x[g[1],j] for j in JJ])) < t[g[2]]
                    t[g[2]] = t[g[1]] - D[g[1]]*(1 - sum([ee[j]*x[g[1],j] for j in JJ]));
                end
            else
                if t[g[1]] - D[g[1]] < t[g[2]]
                    t[g[2]] = t[g[1]] - D[g[1]];
                end
            end
        end
    end
    tv = zeros(length(II));
    for i in II
        tv[i] = -t[i];
    end
    return tv;
end
