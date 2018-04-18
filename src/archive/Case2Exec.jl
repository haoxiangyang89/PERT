# this is the collection of executions for case 2: the disruption affects the remaining events
#cd("C:\\Documents\\Git\\PERT\\");
include("Case2Func_Cplex.jl");
InputAdd = "test_Input_graph_Full.csv";
D,r,H,b,B,ee,II,JJ,M,SS,GG,dH,dR,p = readIn(InputAdd);

#mf = fullExt(D,r,H,b,B,ee,II,JJ,SS,GG,dH,dR,p,M);
#solve(mf);
#mfObj = getobjectivevalue(mf);
#mfX = getvalue(mf[:x])[:,:,1];
#mft = getvalue(mf[:t])[:,1];

# intialize with the solutions
mi = makeMaster(D,r,H,b,B,ee,II,JJ,SS,GG,dH,dR,p);

noS = length(SS) - 1;
cumS = Dict{Any,Any}();
rscen = Dict{Any,Any}();
for s in SS[2:length(SS)]
    rscen[s] = Dict{Any,Any}();
    for kr in keys(r)
        if kr[2] == s
            rscen[s][kr[1]] = r[kr];
        end
    end
end

# Step 0: initialization
k = 1;
l = zeros(length(SS));
ub = Any[9999999];
lb = Any[0];
keepIter = true;
ϵ = 0.01;
counter = 0;
maxIter = 3;
cuts = Dict{Any,Any}();
xrec = [];
trec = [];
xfinrec = [];
tfinrec = [];
ztight = 0;
solver = CplexSolver(CPXPARAM_Threads = 3, CPXPARAM_ScreenOutput = 0);

while keepIter
    # obtain the current master solutions
    solve(mi);
    miObj = getobjectivevalue(mi);
    miX = getvalue(mi[:x])[:,:];
    push!(xrec,miX);
    mit = getvalue(mi[:t])[:];
    push!(trec,mit);
    cumuF = p*maximum(mit);

    # obtain the lower bound
    if miObj > lb[length(lb)]
        push!(lb,miObj);
    else
        temp = lb[length(lb)];
        push!(lb,temp);
    end

    # generate Lagrangian cuts
    for s in SS[2:length(SS)]
        I1,I2,I3,tau = obtainIs(miX,mit,D,rscen[s],H[s],b,B,ee,II,JJ,GG);
        msI = subInt(miX,mit,D,tau,b,B,ee,II,I1,I2,I3,JJ,GG);
        solve(msI);
        msIObj = getobjectivevalue(msI);
        cumS[s] = msIObj;
        cumuF += (1-p)/noS*msIObj;

        # solve the lagrangian prime to generate the cuts
        ppi,lambda1,lambda2,lambda,zlagp = solveSub(miX,mit,D,rscen[s],H[s],b,B,ee,II,JJ,GG,cumS[s]);
        mlag = subLag(miX,mit,D,rscen[s],H[s],M[s],b,B,ee,II,I1,I2,I3,JJ,GG,lambda1,lambda2,ppi);
        solve(mlag);
        zlag = getobjectivevalue(mlag);
        if abs(zlagp - zlag) < ϵ
            ztight += 1;
        end
      #  xreal1 = getvalue(mlag[:xre1]);
      #  xreal2 = getvalue(mlag[:xre2]);
      #  xreal3 = getvalue(mlag[:xre3]);
      #  ppit,lambda1t,lambda2t,lambdat,zlag = solveSub2(miX,mit,D,rscen[s],H[s],b,B,ee,II,JJ,GG,cumS[s],M);
        l[s] += 1;

        # append Lagrangian cuts to the master program
        cuts[s,l[s]] = @constraint(mi,mi[:θ][s] >= sum{ppi[i]*(mi[:t][i] - mit[i]),i in II;i in keys(ppi)}
            + sum{lambda[i,j]*(mi[:x][i,j] - miX[i,j]), i in II, j in JJ;(i,j) in keys(lambda)} + zlag);
    end

    # obtain the upper bound
    if cumuF < ub[length(ub)]
        push!(ub,cumuF);
        push!(xfinrec,miX);
        push!(tfinrec,mit);
    else
        push!(ub,ub[length(ub)]);
    end

    # decides whether to keep the iteration
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

# obtain the final solution and evaluation
fvfin = ub[length(ub)];
xfin = xfinrec[length(xfinrec)];
tfin = tfinrec[length(tfinrec)];

# print out the result
println("---------------------lb--------------------");
println(lb);
println("---------------------ub--------------------");
println(ub);
println("---------------------x---------------------");
println(xfin);
println("---------------------t---------------------");
println(tfin);
println("---------------------f---------------------");
println(fvfin);
