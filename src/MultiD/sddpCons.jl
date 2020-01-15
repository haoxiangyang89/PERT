# multiple disruption SDDP construction
using SDDP, Gurobi;
using Distributions,HDF5,JLD,DelimitedFiles;
include("header.jl");

# initiate the data
filePath = "/Users/haoxiangyang/Desktop/PERT_tests/current/11";
hMax = 10000;
Ωsize = 1;
pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize,hMax);
nStage = 3;
integrality_handler = SDDP.ContinuousRelaxation();
# integrality_handler = SDDP.SDDiP();
# disData is a dictionary, keys are the stages
disData = Dict();
Ωn = 10;
for n in 1:nStage
    disData[n],Ω = autoUGen(nameH, Hparams, nameD, dparams, Ωn, 1);
end

model = SDDP.LinearPolicyGraph(
        stages = nStage,
        sense = :Min,
        lower_bound = 0.0,
        optimizer = with_optimizer(Gurobi.Optimizer),
        integrality_handler = integrality_handler,
    ) do subproblem, stage
    set_parameter(subproblem, "IntFeasTol", 1e-8)
    set_parameter(subproblem, "FeasibilityTol", 1e-8)
    set_parameter(subproblem, "OutputFlag", 0)
    add_state_variables(subproblem, pData)
    if stage == 1
        # create first stage
        create_first_stage(subproblem, pData)
    else
        # create later stages
        create_n_stage(subproblem, stage, pData, disData[stage],10000)
    end
end

SDDP.train(model; iteration_limit = 30);
