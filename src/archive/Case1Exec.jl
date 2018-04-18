# this is the collection of executions for case 1: the disruption does not affect the events already started

include("Case1Func.jl");
include("def.jl");
# read in data of scenarios and the activity network
# Semi1 is the one with fixed disruption time
# Semi2 is the one with fixed disruption magnitude
# Fixed is the one with both disruption time and disruption magnitude fixed
# InputAdd = "test_Input_graph_Fixed.csv";
#InputAdd = "test_Input_graph_Semi1.csv";
#D,d,H,b,B,ee,II,JJ,M,SS,GG,dH,dR,p = readIn(InputAdd);
jldAdd = "initFile.jld";
D,d,H,b,B,ee,II,JJ,M,SS,GG,p = loadInit(jldAdd);
ϵ = 1e-4;
maxIter = 10;

# prepare the scenario data: each scenario's information with each activity's duration after disruption
noS = length(SS) - 1;
dscen = Dict{Any,Any}();
for s in SS[2:length(SS)]
    dscen[s] = Dict{Any,Any}();
    for kr in keys(d)
        if kr[2] == s
            dscen[s][kr[1]] = d[kr];
        end
    end
end

# step 1: initialization
nodeList = [];
# create the master problem
mp = createMaster(D,d,b,B,ee,II,JJ,SS,GG,p);
# create the initial node
bSetIni = Dict();
bSignSetIni = Dict();
for s in SS[2:length(SS)]
  bSetIni[s] = [];
  bSignSetIni[s] = [];
end
iniNode = nodeType(0,0,Inf,mp,bSetIni,bSignSetIni);
nodeList = [iniNode];
# the universal lower bound/upper bound
totalLB = 0;
totalUB = Inf;
nodeCount = 0;
tSolUB = Dict();
xSolUB = Dict();

while nodeList != []
  # solve the master problem for t/x
  currentNode = nodeList[1];
  if currentNode.lbCost < totalUB
    mpi = currentNode.mp;
    keepIter = true;
    noImprove = 0;
    while keepIter
      pastLB = copy(currentNode.lbCost);
      pastUB = copy(currentNode.ubCost);
      solve(mpi);
      mpObj = getobjectivevalue(mpi);
      # update the node lower bound
      if mpObj > currentNode.lbCost
        currentNode.lbCost = mpObj;
      end

      # collect the solution from the master program
      tSol = Dict();
      xSol = Dict();
      for i in II
        tSol[i] = getvalue(mpi[:t][i]);
        for j in JJ
          xSol[i,j] = getvalue(mpi[:x][i,j]);
        end
      end

      # for each subproblem, solve the subproblem corresponding to the current node
      tempUB = p*(getvalue(mpi[:tN]));
      for s in SS[2:length(SS)]
        I1,I2 = obtainIs(xSol,tSol,H[s],II);

        # obtain the upper bound of the this subproblem
        msI = subInt(xSol,tSol,D,b,B,ee,II,I1,I2,dscen[s],JJ,GG);
        solve(msI);
        msIObj = getobjectivevalue(msI);
        tempUB += (1-p)*msIObj/(length(SS)-1);

        # generate Lagrangian cuts
        #πle,λle,zlagp = solveSub(xSol,tSol,D,dscen[s],M[s],H[s],b,B,ee,II,JJ,GG,msIObj,currentNode.bSet[s],currentNode.bSignSet[s]);
        πle,λle,zlagp = solveSubI(xSol,tSol,D,dscen[s],M[s],H[s],b,B,ee,II,JJ,GG,msIObj,currentNode.bSet[s],currentNode.bSignSet[s]);
        mlag = subLag(xSol,tSol,D,dscen[s],H[s],M[s],b,B,ee,II,I1,I2,JJ,GG,πle,λle);
        for i in 1:length(currentNode.bSet[s])
          if currentNode.bSignSet[s][i] == 1
            @constraint(mlag,mlag[:t][currentNode.bSet[s][i]] <= H[s] - 1e-4);
          elseif currentNode.bSignSet[s][i] == 2
            @constraint(mlag,mlag[:t][currentNode.bSet[s][i]] >= H[s]);
          end
        end
        solve(mlag);
        zlag = getobjectivevalue(mlag);

        # append Lagrangian cuts to the current node
        cutTemp = generateCut(s,πle,λle,zlag,xSol,tSol,I1,JJ);
        mpi = appendLCcuts(mpi,cutTemp,II,JJ);
      end
      # update the upper bound
      if tempUB < currentNode.ubCost
        currentNode.ubCost = tempUB;
        tSolUB = copy(tSol);
        xSolUB = copy(xSol);
      end
      # whether to end generating Lagrangian cuts
      if currentNode.ubCost - currentNode.lbCost >= ϵ*currentNode.lbCost
        if (currentNode.ubCost - currentNode.lbCost) < (pastUB - pastLB)
          noImprove = 0;
        else
          noImprove += 1;
        end
        if noImprove >= maxIter
          keepIter = false;
        end
      else
        keepIter = false;
      end
    end
    # update the universal upper bound
    if currentNode.ubCost < totalUB
      totalUB = currentNode.ubCost;
    end

    # decide which scenario time and which activity to branch on
    solve(mpi);
    tSol = Dict();
    xSol = Dict();
    for i in II
      tSol[i] = getvalue(mpi[:t][i]);
      for j in JJ
        xSol[i,j] = getvalue(mpi[:x][i,j]);
      end
    end
    zDiffMax = 0;
    zDiffMaxIndex = 0;
    for s in SS[2:length(SS)]
      θhat = getvalue(mpi[:θ][s]);
      I1,I2 = obtainIs(xSol,tSol,H[s],II);

      # obtain the upper bound of the this subproblem
      msI = subInt(xSol,tSol,D,b,B,ee,II,I1,I2,dscen[s],JJ,GG);
      solve(msI);
      θbar = getobjectivevalue(msI);
      if zDiffMax < (θbar - θhat) - 1e-5
        zDiffMax = (θbar - θhat);
        zDiffMaxIndex = s;
      end
    end

    if zDiffMax > 0
      # branch into two problems: branch on the last event before the disruption time in scenario zDiffMaxIndex
      # find the last event started before the disruption time
      bestE = +Inf;
      bestEIndex = 0;
      for i in II
        # if we have not branched on this scenario time on this activity
        # if (abs(tSol[i] - H[zDiffMaxIndex]) < bestE)&(!(i in currentNode.bSet[zDiffMaxIndex]))
        if (abs(tSolUB[i] - H[zDiffMaxIndex]) < bestE)&(!(i in currentNode.bSet[zDiffMaxIndex]))
          # bestE = abs(tSol[i] - H[zDiffMaxIndex]);
          bestE = abs(tSolUB[i] - H[zDiffMaxIndex]);
          bestEIndex = i;
        end
      end
      nodeCount += 1;
      bSetTemp = deepcopy(currentNode.bSet);
      bSignSetTemp1 = deepcopy(currentNode.bSignSet);
      bSignSetTemp2 = deepcopy(currentNode.bSignSet);
      push!(bSetTemp[zDiffMaxIndex],bestEIndex);

      push!(bSignSetTemp1[zDiffMaxIndex],1);
      mpBran1 = copy(mpi);
      mpBran1 = appendBNBcuts(mpBran1,bestEIndex,H[zDiffMaxIndex],1);
      node1 = nodeType(nodeCount,currentNode.lbCost,currentNode.ubCost,mpBran1,bSetTemp,bSignSetTemp1);
      push!(nodeList,node1);

      nodeCount += 1;
      mpBran2 = copy(mpi);
      push!(bSignSetTemp2[zDiffMaxIndex],2);
      mpBran2 = appendBNBcuts(mpBran2,bestEIndex,H[zDiffMaxIndex],2);
      node2 = nodeType(nodeCount,currentNode.lbCost,currentNode.ubCost,mpBran2,bSetTemp,bSignSetTemp2);
      push!(nodeList,node2);
    end
  end
  # remove the current node
  shift!(nodeList);
end
