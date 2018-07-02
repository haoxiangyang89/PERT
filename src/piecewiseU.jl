import Base.Random.rand

struct piecewiseUniformSampler <: Sampleable{Univariate,Continuous}
    endPoints::Vector{Float64}
    mass::Vector{Float64}
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
