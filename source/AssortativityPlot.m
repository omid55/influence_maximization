% Omid55
%% Initiallization
load('AssortativeNetsData');
N = size(nets{1},1);
InformedAgentsSize = ceil(N/10);
MaximumSimulationSteps = 1000;
mu = 0.9;
beta = 0.1;
Randomness = 0.01;
alpha = 0;
runs = 100;

%% Simulation
% Betweenness & Closeness Calculation
list = 0.9:-0.2:-0.9;
finalOp = zeros(runs,length(list));
for num = 1 : length(list)
    fprintf(1,'Algorithm Num %d\n',num);
    A = nets{num};
    sp = sparse(A);
    betweennesses = betweenness_centrality(sp);
    [D ~] = all_shortest_paths(sp,struct('algname','floyd_warshall'));
    closenesses = (N-1) ./ sum(D);
    clear('sp');

    for run = 1 : runs
        x = -1 + 2 * rand(N,1);
        [meanMajOpininos4,d4,fd4,fol4] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,4,Randomness,alpha,beta,betweennesses,closenesses);
        finalOp(run,num) = meanMajOpininos4(end);
    end
end

errorbar(list,mean(finalOp),2*std(finalOp)/sqrt(runs),'LineWidth',2);
save('AssortativeFinalData');
    
