# this is the collection of executions for case 1: the disruption does not affect the events already started

include("Case1Func.jl");
# read in data of scenarios and the activity network
InputAdd = "test_Input_graph_Full.csv";
D,d,H,b,B,ee,II,JJ,M,SS,GG,dH,dR,p = readIn(InputAdd);

# step 1: initialization
nodeList = [];
# create the initial node
currentNode = 

# create the master problem
mp = createMaster(D,d,H,b,B,ee,II,JJ,SS,GG,dH,dR,p);
# solve the master problem for t/x

# solve the subproblem corresponding to the current node

# generate Lagrangian cuts

# append Lagrangian cuts to the current node

# branch into two problems

# remove the current node
