%Omid55
%Proposed Algorithm run
function [ cost ] = CostFunction( chromosome,x,N,InformedAgentsSize,MaximumSimulationSteps,XD,UI,alpha,net,EPS )

% decoding
mu = chromosome(1);
muPrime = chromosome(2);
c1 = chromosome(3);
c2 = chromosome(4);

cost = 0;
runs = 5;
for run=1:runs
    net1 = containers.Map(net.keys,net.values);
    [meanAllOpinions,meanAgOpininos,meanMajOpininos] = SimulationMethod(x,N,InformedAgentsSize,MaximumSimulationSteps,mu,muPrime,XD,UI,alpha,net1,3,c1,c2,EPS);
    cost = cost + 1 - meanMajOpininos(end);
end
cost = cost / runs;

end

