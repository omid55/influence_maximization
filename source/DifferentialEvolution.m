%Omid55
%DE
function [ bestChromosomeCost,bestChromosome ] = DifferentialEvolution( x,N,InformedAgentsSize,MaximumSimulationSteps,XD,UI,alpha,net,EPS )

disp('DE');
dimension = 4;
MaxIteration = 10;
PopulationSize = 5;
Beta = 0.99;
Pr = 0.99;
Nv = 1;   % number of difference vectors    
de_it = 0;
de_population = rand(PopulationSize,dimension);
de_population(:,2) = 0.7 + de_population(:,2)*0.3;
de_population(:,4) = 0.7 + de_population(:,4)*0.3;

%%%
costs = [];
%%%

% % % costs = zeros(PopulationSize,1);
% % % for i=1:size(de_population,1)
% % %     costs(i) = CostFunction(de_population(i,:),x,N,InformedAgentsSize,MaximumSimulationSteps,XD,UI,alpha,net,sp,EPS);
% % % end

de_bests = [];
de_means = [];
while de_it < MaxIteration
    lastDe = de_population;

    de_it = de_it + 1;
    fprintf(1,'Differential Evolution Iteration %d\n',de_it);
% % %     de_bests = [de_bests; min(costs)];
% % %     de_means = [de_means; mean(costs)];
% % % 
% % %     plot(1:de_it,de_bests,1:de_it,de_means);
% % %     legend('Min of Costs','Mean of Costs');
% % %     xlabel('Generations');
% % %     ylabel('Fitnesses');
% % %     title('Differential Evolution');

    [de_population] = CreateTrialVectorAndCrossOver(de_population,costs,Beta,Pr,Nv,x,N,InformedAgentsSize,MaximumSimulationSteps,XD,UI,alpha,net,EPS);

% % %     for i=1:size(de_population,1)
% % %         costs(i) = CostFunction(de_population(i,:),x,N,InformedAgentsSize,MaximumSimulationSteps,XD,UI,alpha,net,sp,EPS);
% % %     end

    if norm(de_population - lastDe) < EPS
        break;
    end
end

%%%
for i=1:size(de_population,1)
    costs(i) = CostFunction(de_population(i,:),x,N,InformedAgentsSize,MaximumSimulationSteps,XD,UI,alpha,net,EPS);
end
%%%

bestIdx = find(costs == min(costs));
bestChromosomeCost = costs(bestIdx);
bestChromosome = de_population(bestIdx,:);
    
end

