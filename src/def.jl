# define the customized variable format

type nodeType
  nodeID :: Int64
  lbCost :: Float64
  ubCost :: Float64
  mp :: JuMP.Model
  bSet :: Dict{Any,Any}
  bSignSet :: Dict{Any,Any}
end

type LagCut
  s :: Int64
  # coefficients for t
  πl :: Dict{Any,Any}
  # coefficients for t
  λl :: Dict{Any,Any}
  # constant term
  z :: Float64
end
