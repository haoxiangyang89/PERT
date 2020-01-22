import Base.rand,Statistics.mean

struct piecewiseUniformSampler <: Sampleable{Univariate,Continuous}
    endPoints::Array{Float64,1}
    mass::Array{Float64,1}
end

struct CategoricalSamplerNew <: Sampleable{Univariate,Discrete}
    mass::Array{Float64,1}
    category::Array{Any,1}
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

function rand(s::CategoricalSamplerNew)
    s.category[rand(Categorical(s.mass))]
end

function mean(s::CategoricalSamplerNew)
    # calculate the meann of the piecewise sampler
    y = 0
    for i in 1:length(s.mass)
        y += s.mass[i]*s.category[i];
    end
    y
end
