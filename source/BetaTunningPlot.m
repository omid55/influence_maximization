%% Different Beta
close all;
clc;
clear;

%% Initiallization
N = 1000;
MaximumSimulationSteps = 350;
mu = 0.9;
runs = 20;
Randomness = 0.01;
averageDegree = 12;
A = BAnet(N,averageDegree/2,averageDegree/2);
sp = sparse(A);
InformedAgentsSize = ceil(N / 10);
alpha = 0;

%% Simulation
list = 0:0.1:1;
leng = MaximumSimulationSteps + 1;
meanMajOpininos1Avg = zeros(runs,length(list));
x = -1 + 2 * rand(N,1);
for run=1:runs
    cnt = 1;
    fprintf(1,'Algorithm Run %d\n',run);
    for beta = list
        [meanMajOpininos1,~,~,~] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,1,Randomness,alpha,beta);
        meanMajOpininos1Avg(run,cnt) = meanMajOpininos1(end);
        cnt = cnt + 1;
    end
end

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 6 4]);
h = errorbar(list,mean(meanMajOpininos1Avg),std(meanMajOpininos1Avg)/sqrt(runs),'LineWidth',2);
xlabel('\bf\mu');
ylabel('\bfOpinion of population at The End of Simulation');
