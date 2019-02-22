import Base.Random.rand,Base.mean

struct piecewiseUniformSampler <: Sampleable{Univariate,Continuous}
    endPoints::Vector{Float64}
    mass::Vector{Float64}
end

struct CategoricalSampler <: Sampleable{Univariate,Discrete}
    mass::Vector{Float64}
    category::Vector{Float64}
end

function rand(s::piecewiseUniformSampler)
    l = rand();
    pieces = length(s.mass);
    remain = l;
    j = 0;
    while remain >= 0
        j += 1;
        remain -= s.mass[j];
    end
    x = (remain + s.mass[j])/s.mass[j]*(s.endPoints[j+1] - s.endPoints[j]) + s.endPoints[j];
    x
end

function mean(s::piecewiseUniformSampler)
    # calculate the meann of the piecewise sampler
    y = 0
    for i in 1:length(s.mass)
        y += s.mass[i]*(s.endPoints[i+1] + s.endPoints[i])/2;
    end
    y
end

function rand(s::CategoricalSampler)
    s.category[rand(Categorical(s.mass))]
end
