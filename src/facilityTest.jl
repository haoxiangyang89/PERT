function masterFacility(f,c,d,M,I,J)
    mp = Model(solver = GurobiSolver());
    @variable(mp, x[i in I], Bin);
    @variable(mp, θ >= 0);
    @objective(mp, Min, sum(f[i]*x[i] for i in I) + θ);
    return mp;
end

function subFacility(f,c,cz,d,M,I,J,xhat)
    sp = Model(solver = GurobiSolver());
    @variable(sp, 0 <= xs[i in I] <= 1);
    @variable(sp, y[i in I,j in J] >= 0);
    @variable(sp, z[j in J] >= 0);
    @constraint(sp, supply[i in I], sum(y[i,j] for j in J) <= M[i]*xs[i]);
    @constraint(sp, demand[j in J], sum(y[i,j] for i in I) + z[j] >= d[j]);
    @constraint(sp, xeqn[i in I], xs[i] == xhat[i]);
    @objective(sp, Min, sum(sum(c[i,j]*y[i,j] for j in J) for i in I) + sum(cz[j]*z[j] for j in J));

    solve(sp);
    πDict = Dict();
    for i in I
        πDict[i] = getdual(sp[:xeqn][i]);
    end
    v = getobjectivevalue(sp);
    return πDict,v;
end

function addFacCuts(mp,πDict,v,xhat,I)
    @constraint(mp, mp[:θ] >= v + sum(πDict[i]*(mp[:x][i] - xhat[i]) for i in I));
    return mp;
end

function subFacility_Lag(f,c,cz,d,M,I,J,xhat,πCurr)
    sp = Model(solver = GurobiSolver());
    @variable(sp, 0 <= xs[i in I] <= 1);
    @variable(sp, y[i in I,j in J] >= 0);
    @variable(sp, z[j in J] >= 0);
    @constraint(sp, supply[i in I], sum(y[i,j] for j in J) <= M[i]*xs[i]);
    @constraint(sp, demand[j in J], sum(y[i,j] for i in I) + z[j] >= d[j]);
    @objective(sp, Min, sum(sum(c[i,j]*y[i,j] for j in J) for i in I) + sum(cz[j]*z[j] for j in J) + sum(πCurr[i]*(xs[i] - xhat[i]) for i in I));

    solve(sp);
    v = getobjectivevalue(sp);
    xsDict = Dict();
    for i in I
        xsDict[i] = getvalue(sp[:xs][i]);
    end
    return v,xsDict;
end

function LagrangianFac(f,c,cz,d,M,I,J,xhat)
    sp = Model(solver = GurobiSolver());
    @variable(sp, 0 <= xs[i in I] <= 1);
    @variable(sp, y[i in I,j in J] >= 0);
    @variable(sp, z[j in J] >= 0);
    @constraint(sp, supply[i in I], sum(y[i,j] for j in J) <= M[i]*xs[i]);
    @constraint(sp, demand[j in J], sum(y[i,j] for i in I) + z[j] >= d[j]);
    @constraint(sp, xeqn[i in I], xs[i] == xhat[i]);
    @objective(sp, Min, sum(sum(c[i,j]*y[i,j] for j in J) for i in I) + sum(cz[j]*z[j] for j in J));

    solve(sp);
    πDict = Dict();
    for i in I
        πDict[i] = getdual(sp[:xeqn][i]);
    end
    # start from the Benders cut dual
    for k in 1:10
        v,xsDict = subFacility_Lag(f,c,cz,d,M,I,J,xhat,πDict);
        for i in I
            πDict[i] = πDict[i] + (3000 - v)/sum((xsDict[i] - xhat[i])^2 for i in I)*(xsDict[i] - xhat[i]);
        end
    end
    v,xsDict = subFacility_Lag(f,c,cz,d,M,I,J,xhat,πDict);
    return v,πDict;
end

f = [300,300];
M = [300,300];
c = [1 2 10;10 2 1];
d = [100,100,100];
I = 1:2;
J = 1:3;
# Benders iterations
mp = masterFacility(f,c,d,M,I,J);
solve(mp);
xhat = Dict();
for i in I
    xhat[i] = getvalue(mp[:x][i]);
end
πDict,v = subFacility(f,c,cz,d,M,I,J,xhat);
mp = addFacCuts(mp,πDict,v,xhat,I);


# test the optimal cutting plane generation
mpt = Model(solver = GurobiSolver());
@variable(mpt, x[i in I], Bin);
@variable(mpt, θ >= 0);
@objective(mpt, Min, sum(f[i]*x[i] for i in I) + θ);
@constraint(mpt, θ >= 3000 - 1700*x[1] - 1700*x[2]);
#@constraint(mpt, θ >= 2200 - 900*x[1] - 900*x[2]);
solve(mpt,relaxation = true);
