# this is a test file of second stage for PERT using a pilot example

using JuMP, Cbc;

m = Model();
I = 1:6;
@variable(m,0 <= x1[1:2] <= 1, Bin);
@variable(m,x2[3:6],Bin);
@variable(m,t[I] >= 0);

@constraint(m,tbinding[i in 1:2],t[i] == 0);
@constraint(m,t31,t[3]-t[1]>=20*(1-0.5*x1[1]));
@constraint(m,t32,t[3]-1.5*t[2]>=45*(1-0.5*x1[2])-5);
@constraint(m,t43,t[4]-t[3]>=10.5*(1-0.5*x2[3]));
@constraint(m,t53,t[5]-t[3]>=10.5*(1-0.5*x2[3]));
@constraint(m,t64,t[6]-t[4]>=18*(1-0.5*x2[4]));
@constraint(m,t64,t[6]-t[5]>=12*(1-0.5*x2[5]));
@constraint(m,budget,sum{x1[i],i in 1:2} + sum{x2[i],i in 3:6}<=2);

solver = CbcSolver();
valueMat = zeros(41,41);
optV = -99999;
lambdO1 = 0;
lambdO2 = 0;
for lambda1 in 0:20
    for lambda2 in -20:20
        @objective(m,Min,t[6]-lambda1*(x1[1]-1)-lambda2*(x1[2]-1));
        solve(m);
        currentOpt=getobjectivevalue(m);
        #valueMat[lambda1+21,lambda2+21]=currentOpt;
        if optV < currentOpt
            optV=currentOpt;
            lambdO1=lambda1;
            lambdO2=lambda2;
        end
    end
end

lambda1 = 4.1;
lambda2 = -5.1;
m1 = Model();
@variable(m1,0 <= x1[1:2] <= 1);
@variable(m1,x2[3:6],Bin);
@variable(m1,t[1:6] >= 0);

@constraint(m1,t31,t[3]-t[1]>=20*(1-0.5*x1[1]));
@constraint(m1,t32,t[3]-1.5*t[2]>=45*(1-0.5*x1[2])-5);
@constraint(m1,t43,t[4]-t[3]>=10.5*(1-0.5*x2[3]));
@constraint(m1,t53,t[5]-t[3]>=10.5*(1-0.5*x2[3]));
@constraint(m1,t64,t[6]-t[4]>=18*(1-0.5*x2[4]));
@constraint(m1,t64,t[6]-t[5]>=12*(1-0.5*x2[5]));
@constraint(m1,budget,sum{x1[i],i in 1:2} + sum{x2[i],i in 3:6}<=2);

@objective(m1,Min,t[6]-lambda1*(x1[1]-1)-lambda2*(x1[2]-0));
@constraint(m1,tbinding[i in 1:2],t[i] == 0);

fvalue = zeros(30,30);
for i in 0:29
    for j in 0:29
        JuMP.setRHS(tbinding[1],i);
        JuMP.setRHS(tbinding[2],j);
        solve(m1);
        fvalue[i+1,j+1] = getobjectivevalue(m1);
    end
end
