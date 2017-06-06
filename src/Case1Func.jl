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
        M[s] = sum(values(D))*maximum([r[i,s] for i in II]);
    end
    H[1] = maximum(values(M));
    M[1] = 2*sum(values(D));
    # set of scenarios
    SS = 1:S;

    return D,d,H,b,B,ee,II,J,M,SS,G,distrH,disRList,p
end
