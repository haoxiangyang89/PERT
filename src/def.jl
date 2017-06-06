# define the customized variable format

type nodeType
  nodeID :: Int64
  lbCost :: Float64
  ubCost :: Float64
  # Lagrangian cuts
  LC :: Array{Any,1}
  # branch-and-bound constraints
  SC :: Array{Any,1}
end
