# This is the parameter definition code
# define the customized variable format

# define the project information
type pInfo
  II :: Array{Any,1}
  Ji :: Dict{Any,Any}
  D :: Dict{Any,Any}
  b :: Dict{Any,Any}
  eff :: Dict{Any,Any}
  B :: Float64
  p0 :: Float64

  K :: Array{Any,1}
  Pre :: Dict{Any,Any}
  Succ :: Dict{Any,Any}
end

# define the disruption information
type disInfo
  H :: Float64
  d :: Dict{Any,Any}
  prDis :: Float64
end

# define the node type within the B&C tree
type nodeType
  lbCost :: Float64
  mp :: JuMP.Model
  brInfo :: Array{Any,2}
end

type nodeTypeP
  lbCost :: Float64
  brInfo :: Array{Any,2}
  tmaxD :: Dict{Any,Any}
  tminD :: Dict{Any,Any}
  cutSet :: Array{Any,1}
  state :: Bool
end

# define the cut type for B&C
type cutType
  # coefficients for t
  π :: Dict{Any,Any}
  # coefficients for t
  λ :: Dict{Any,Any}
  # constant term
  v :: Float64
end

# define the partition type
type partType
  startH :: Int64
  endH :: Int64
end
