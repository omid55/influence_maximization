%% Different Degree Distribution Ploting
close all;
clc;
clear;

%% Initiallization
N = 100;
MaximumSimulationSteps = 1000;
mu = 0.9;
runs = 100;
Randomness = 0.01;
scale = 0.7;
alpha = 0;
averageDegree = 12;
A = BAnet(N,averageDegree/2,averageDegree/2);
sp = sparse(A);
beta = 0.1;


%% Simulation
leng = MaximumSimulationSteps + 1;
val2Lev = zeros(runs,20);
val1Lev = zeros(runs,20);

for run=1:runs
    x = -1 + 2 * rand(N,1);
    cnt = 1;
    fprintf(1,'Algorithm Run %d\n',run);
    for coef = 0.01:0.01:0.2
        InformedAgentsSize = floor(N*coef);
        [meanMajOpininos1Lev,~,~,~] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,1,Randomness,alpha,beta);
        [meanMajOpininos4Lev,~,~,~] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,4,Randomness,alpha,beta);
        val2Lev(run,cnt) = meanMajOpininos1Lev(end);
        val1Lev(run,cnt) = meanMajOpininos4Lev(end);
        cnt = cnt + 1;
    end
end

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 6 4]);
h = errorbar(0.01:0.01:0.2,mean(val1Lev),std(val1Lev)/sqrt(runs));
set(h,'Color','b');
hold on;
h = errorbar(0.01:0.01:0.2,mean(val2Lev),std(val2Lev)/sqrt(runs));
set(h,'Color','r');
leg = legend('1 Level','2 Level');
set(leg,'FontSize',7);
xlabel('\bf\mu');
ylabel('\bfOpinion of population at The End of Simulation');

save('InformedAgentsRatio');
