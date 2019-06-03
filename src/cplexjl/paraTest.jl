@everywhere function solveModel(data)
    mp = Model(solver = GurobiSolver(Threads = data[4]));
    a = data[1];
    m,n = size(a);
    b = data[2];
    c = data[3];
    @variable(mp,x[1:n] >= 0,Int);
    @constraint(mp, a*x .<= b);
    @objective(mp, Max, sum(c.*x));
    solve(mp);
    return getobjectivevalue(mp);
end

function f_pmap(f, lst)
    npList = procs()            # Number of processes available.
    n  = length(lst)
    results = Vector{Any}(n) # Where we will write the results. As we do not know
                             # the type (Integer, Tuple...) we write "Any"
    i = 1
    nextidx() = (idx = i; i += 1; idx) # Function to know which is the next work item.
                                       # In this case it is just an index.
    @sync begin # See below the discussion about all this part.
        for p in npList
            if p != myid() || length(npList) == 1
                @async begin
                    while true
                        idx = nextidx()
                        println(p);
                        if idx > n
                            break
                        end
                        results[idx] = remotecall_fetch(f, p, lst[idx])
                        println(results[idx]);
                    end
                end
            end
        end
    end
    results
end

function pTest()
    np = nprocs()            # Number of processes available.
    nodeList = ["1"];
    @sync begin
       while nodeList != []
           for p = 1:np
               if p != myid() || np == 1
                   @async begin
                       println(p);
                       while true
                           if length(nodeList) == 0
                               break
                           end
                           b = pop!(nodeList);
                           remotecall_fetch(println,b);
                           if rand() >= 0.5
                               push!(nodeList,"$(b)1");
                               push!(nodeList,"$(b)2");
                           end
                       end
                   end
               end
           end
       end
    end
end
