InputAdd1 = "eval_Full.csv";
D,r1,H1,b,B,ee,II,JJ,M1,SS1,GG,dH1,dR1,p1 = readIn(InputAdd1);
noS1 = length(SS1) - 1;
cumS1 = Dict{Any,Any}();
rscen1 = Dict{Any,Any}();
for s in SS1[2:length(SS1)]
    rscen1[s] = Dict{Any,Any}();
    for kr in keys(r1)
        if kr[2] == s
            rscen1[s][kr[1]] = r1[kr];
        end
    end
end
cumuF1 = p1*maximum(tfin);
for s in SS1[2:length(SS1)]
    I1,I2,I3,tau = obtainIs(xfin,tfin,D,rscen1[s],H1[s],b,B,ee,II,JJ,GG);
    msI1 = subInt(xfin,tfin,D,tau,b,B,ee,II,I1,I2,I3,JJ,GG);
    solve(msI1);
    cumS1[s] = getobjectivevalue(msI1);
    cumuF1 += getobjectivevalue(msI1)*(1-p1)/noS1;
end
