# this is the illustration example Shabbir used in ICSP 2016

using JuMP,Cbc;

m = Model();
solver = CbcSolver();

@variable(m,0<=y1<=2,Int);
@variable(m,0<=y2<=3);
@variable(m,0<=z<=1);

@constraint(m,2*y1+y2>=3*z);
L = -10:10;
for l in L
    @objective(m,Min,y1+y2-(z-0)*l);
    solve(m);
    println(l);
    println(getobjectivevalue(m));
#    print("The optimal value for lambda = ",l," is "getobjectivevalue(m));
end
