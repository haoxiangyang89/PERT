m = Model(solver = GurobiSolver());
that = 0.8;
#@variable(m, 0 <= that <= 1);
@variable(m,x);
@variable(m,0 <= G <= 1);
@variable(m,0 <= F <= 1.1);
@constraint(m, cons1, 0.5 - 0.5G <= that);
@constraint(m, cons2, 1.1 - 0.6G >= that);
@constraint(m, cons3, x + 4*F -2.2*G + 1.2 == 2*that);
@constraint(m, cons4, F <= that);
@constraint(m, cons5, F <= 1.1*G);
@constraint(m, cons6, F >= that - 1.1*(1 - G));

@objective(m,Min,x);
solve(m);

λ1 = getdual(cons1);
λ2 = getdual(cons2);
λ3 = getdual(cons3);
λ4 = getdual(cons4);
λ5 = getdual(cons5);
λ6 = getdual(cons6);
println(λ1 + λ2 + 2*λ3 + λ4 + λ6);


m = Model(solver = GurobiSolver());
that = 0.4;
#@variable(m, 0 <= that <= 1);
@variable(m,x);
@variable(m,0 <= G <= 1);
@variable(m,0 <= F <= 1.1);
#@constraint(m, cons1, 0.5 - 0.5G <= that);
#@constraint(m, cons2, 1.1 - 0.6G >= that);
@constraint(m, cons3, 4*F -2.2*G + 1 - x == 2*that);
@constraint(m, cons4, F <= that);
@constraint(m, cons5, F <= 1.1*G);
@constraint(m, cons6, F >= that - 0.5*(1 - G));
@constraint(m, cons7, F >= 0.5*G);

@objective(m,Min,x);
solve(m);
λ1 = getdual(cons1);
λ2 = getdual(cons2);
λ3 = getdual(cons3);
λ4 = getdual(cons4);
λ5 = getdual(cons5);
λ6 = getdual(cons6);
println(λ1 + λ2 + 2*λ3 + λ4 + λ6);
