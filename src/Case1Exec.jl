# this is the collection of executions for case 1: the disruption does not affect the events already started

include("Case1Func.jl");
# read in data of scenarios and the activity network
InputAdd = "test_Input_graph_Full.csv";
D,d,H,b,B,ee,II,JJ,M,SS,GG,dH,dR,p = readIn(InputAdd);

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
iniNode = nodeType(0,0,Inf,mp);
nodeList = [iniNode];
# the universal lower bound/upper bound
totalLB = 0;
totalUB = Inf;
nodeCount = 0;

while nodeList != []
  # solve the master problem for t/x
  currentNode = nodeList[1];
  if currentNode.lbCost < totalUB
    solve(currentNode.mp);
    mpObj = getobjectivevalue(currentNode.mp);

    # update the node lower bound and the universal lower bound
    if mpObj > currentNode.lbCost
      currentNode.lbCost = mpObj;
    end
    if mpObj > totalLB
      totalLB = mpObj;
    end

    # collect the solution from the master program
    tSol = Dict();
    xSol = Dict();
    for i in II
      tSol[i] = getvalue(mp.varDict[:t][i]);
      for j in JJ
        xSol[i,j] = getvalue(mp.varDict[:x][i,j]);
      end
    end

    # for each subproblem, solve the subproblem corresponding to the current node
    tempUB = p*mpObj;
    zDiff = Dict();
    zDiffMax = -Inf;
    zDiffMaxIndex = 1;
    branchBool = false;
    mpBran = copy(currentNode.mp);
    for s in SS[2:length(SS)]
      I1,I2 = obtainIs(xSol,tSol,H[s],II);

      # obtain the upper bound of the this subproblem
      msI = subInt(xSol,tSol,D,b,B,ee,II,I1,I2,dscen[s],JJ,GG);
      solve(msI);
      msIObj = getobjectivevalue(msI);
      tempUB += (1-p)*msIObj/(length(SS)-1);

      # generate Lagrangian cuts
      πle,λle,zlagp = solveSub(xSol,tSol,D,dscen[s],H[s],b,B,ee,II,JJ,GG,msIObj);
      mlag = subLag(xSol,tSol,D,dscen[s],H[s],M,b,B,ee,II,I1,I2,I3,JJ,GG,πle,λle);
      solve(mlag);
      zlag = getobjectivevalue(mlag);
      zDiff[s] = zlagp - zlag;
      if zDiff[s] > 1e-4
        branchBool = true;
      end
      if zDiff[s] > zDiffMax
        zDiffMax = zDiff[s];
        zDiffMaxIndex = s;
      end

      # append Lagrangian cuts to the current node
      cutTemp = generateCut(s,πle,λle,zlag,xSol,tSol,II,JJ);
      mpBran = appendLCcuts(mpBran,cutTemp,II,JJ);
    end
    # update the upper bound
    if tempUB < currentNode.ubCost
      currentNode.ubCost = tempUB;
    end
    if tempUB < totalUB
      totalUB = tempUB;
    end

    if branchBool
      # branch into two problems: branch on the last event before the disruption time in scenario zDiffMaxIndex
      # find the last event started before the disruption time
      lastE = -Inf;
      lastEIndex = 0;
      for i in II
        if (tSol[i] < H[zDiffMaxIndex])&(tSol[i] > lastE)
          lastE = tSol[i];
          lastEIndex = i;
        end
      end
      nodeCount += 1;
      mpBran1 = appendBNBcuts(mpBran,lastEIndex,H[zDiffMaxIndex],1);
      node1 = nodeType(nodeCount,currentNode.lbCost,currentNode.ubCost,mpBran1);
      push!(nodeList,node1);
      nodeCount += 1;
      mpBran2 = appendBNBcuts(mpBran,lastEIndex,H[zDiffMaxIndex],2);
      node2 = nodeType(nodeCount,currentNode.lbCost,currentNode.ubCost,mpBran2);
      push!(nodeList,node2);
    end
  end
  # remove the current node
  shift!(nodeList);
end
