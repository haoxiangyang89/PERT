# multiple disruption SDDP construction
using SDDP, Gurobi;
using Distributions,HDF5,JLD,DelimitedFiles;
include("header.jl");

# initiate the data
filePath = "/Users/haoxiangyan/Desktop/PERT_tests/current/14";
hMax = 100000;
Ωsize = 1;
nStage = 3;

pData,disDataSet,nameD,nameH,dparams,Hparams = genData(filePath,Ωsize,hMax);
include("extForm.jl");

disDataOri,Ωlocal = autoUGen(nameH, Hparams, nameD, dparams, 3, 1);
save("disDataTest.jld","disData",disDataOri);
disDataOriRaw = load("disDataTest.jld");
disDataOri = disDataOriRaw["disData"];
global Ω = 1:3;
global HList = [10,40,70];
global dList = [disDataOri[ω].d for ω in Ω];
global prList = [1/3, 1/3, 1/3];

# disData is a dictionary, keys are the stages
disData = Dict();
Ωn = 3;
for n in 1:nStage
    #disData[n],Ω = autoUGen(nameH, Hparams, nameD, dparams, Ωn, 1);
    disData[n] = [];
    for ω in 1:Ωn
        dparams = Dict();
        for i in pData.II
            if i != 0
                dparams[i] = disDataOri[ω].d[i];
            end
        end
        push!(disData[n],disInfo(HList[ω],dparams,prList[ω]))
    end
end

integrality_handler = SDDP.ContinuousRelaxation();
# integrality_handler = SDDP.SDDiP();
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
        create_n_stage(subproblem, stage, pData, disData[stage], 100000)
    end
end

SDDP.train(model; iteration_limit = 300);
lpLb = SDDP.calculate_bound(model);

global mExt = Model(with_optimizer(Gurobi.Optimizer));
mExt = extForm(0,[[],[]],pData,0,Dict(),1,0,0);
optimize!(mExt);
trueObj = objective_value(mExt);

global mExt1 = Model(with_optimizer(Gurobi.Optimizer));
mExt1 = extForm_lin(0,[[],[]],pData,0,Dict(),1,0,0);
optimize!(mExt1);
lpObj = objective_value(mExt1);
