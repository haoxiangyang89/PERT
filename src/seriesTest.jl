# test case with two activities in series
@everywhere using JuMP,Gurobi,CPLEX,Ipopt;
@everywhere using Distributions,HDF5,JLD;

@everywhere include("def.jl");
@everywhere include("readIn.jl");
@everywhere include("master.jl");
@everywhere include("sub.jl");
@everywhere include("cuts.jl");
@everywhere include("detForm.jl");
@everywhere include("extForm.jl");
@everywhere include("ubCalFunc.jl");

x1List = linspace(0,1,11);
t1List = linspace(0,1,11);
t2List = linspace(0,1,11);
vstar = Inf;
solutionSet = [];
objList = [];
for x1s in x1List
    for t1s in t1List
        for t2s in t2List
            if t2s - t1s >= 1.4 - 0.7*x1s
                m = Model(solver = IpoptSolver(linear_solver = "ma27"));
                @variable(m, 0 <= x1 <= 1, start = x1s);
                @variable(m, 0 <= t1 <= 1, start = t1s);
                @variable(m, 0 <= t2 <= 1, start = t2s);
                @constraint(m, t2 - t1 >= 1.4 - 0.7*x1);
                @NLobjective(m, Min, -t1^2/2 -t2^2 + t1*t2 - 0.8*x1*t2 - 0.2*x1*t1 - 0.2326*t1 + 2.3*t2 + x1);
                solve(m);
                push!(objList, getobjectivevalue(m) + 1);
                if getobjectivevalue(m) + 1 < vstar
                    vstar = getobjectivevalue(m) + 1;
                    solutionSet = [getvalue(m[:t1]),getvalue(m[:t2]),getvalue(m[:x1])];
                end
            end
        end
    end
end

vstar2 = Inf;
solutionSet2 = [];
objList2 = [];
for x1s in x1List
    for t1s in t1List
        if t1s + 1.4 - 0.7*x1s >= 1
            m = Model(solver = IpoptSolver(linear_solver = "ma27"));
            @variable(m, 0 <= x1 <= 1, start = x1s);
            @variable(m, 0 <= t1 <= 1, start = t1s);
            @constraint(m, t1 + 1.4 - 0.7*x1 >= 1);
            @NLobjective(m, Min, -t1^2/2 - 0.2*x1*t1 + 0.7674*t1+ 0.2*x1);
            solve(m);
            push!(objList2, getobjectivevalue(m) + 2.3);
            if getobjectivevalue(m) + 2.3 < vstar2
                vstar2 = getobjectivevalue(m) + 2.3;
                solutionSet2 = [getvalue(m[:t1]),getvalue(m[:x1])];
            end
        end
    end
end

# construct SAA test
pInputAdd = "test_2_P.csv";
kInputAdd = "test_2_K.csv";
ϕInputAdd = "test_2_Phi_full.csv";

pData = readInP(pInputAdd,kInputAdd);
nameD,dparams = readInUnc(ϕInputAdd);
disData,Ω = autoUGen("Uniform",[0,1],nameD,dparams,5000,1 - pData.p0);
disData = orderdisData(disData,Ω);
text,xext,fext,gext,mp = extForm_cheat(pData,disData,Ω);
