%% Different Social Powers
close all;
clc;
clear;

%% Initiallization
N = 100;
% MaximumSimulationSteps = 100000;
% mu = 0.9;
MaximumSimulationSteps = 1000;
mu = 0.4;
runs = 10;
Randomness = 0.01;
averageDegree = 12;
% A = BAnet(N,averageDegree/2,averageDegree/2);
% sp = sparse(A);clearclc
sp = ForestFireCreator(N,0.2,0.2);
A = full(sp);
InformedAgentsSize = ceil(N / 10);
beta = 0.1;


%% Simulation
list = 0:0.1:2;
leng = MaximumSimulationSteps + 1;
meanMajOpininos1Avg = zeros(runs,length(list));
meanMajOpininos2Avg = zeros(runs,length(list));
meanMajOpininos3Avg = zeros(runs,length(list));
meanMajOpininos4Avg = zeros(runs,length(list));
meanMajOpininos5Avg = zeros(runs,length(list));
meanMajOpininos6Avg = zeros(runs,length(list));

%% Betweenness & Closeness Calculation
betweennesses = betweenness_centrality(sp);
[D ~] = all_shortest_paths(sp,struct('algname','floyd_warshall'));
closenesses = (N-1) ./ sum(D);

for run=1:runs
    x = -1 + 2 * rand(N,1);
    cnt = 1;
    fprintf(1,'Algorithm Run %d\n',run);
    for alpha = list
        [meanMajOpininos1,~,~,~] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,1,Randomness,alpha,beta,betweennesses,closenesses);    
        [meanMajOpininos2,~,~,~] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,2,Randomness,alpha,beta,betweennesses,closenesses);    
        [meanMajOpininos3,~,~,~] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,3,Randomness,alpha,beta,betweennesses,closenesses);       
        [meanMajOpininos4,~,~,~] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,4,Randomness,alpha,beta,betweennesses,closenesses);    
        [meanMajOpininos5,~,~,~] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,5,Randomness,alpha,beta,betweennesses,closenesses);    
        [meanMajOpininos6,~,~,~] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,6,Randomness,alpha,beta,betweennesses,closenesses);    
        meanMajOpininos1Avg(run,cnt) = meanMajOpininos1(end);
        meanMajOpininos2Avg(run,cnt) = meanMajOpininos2(end);
        meanMajOpininos3Avg(run,cnt) = meanMajOpininos3(end);
        meanMajOpininos4Avg(run,cnt) = meanMajOpininos4(end);
        meanMajOpininos5Avg(run,cnt) = meanMajOpininos5(end);
        meanMajOpininos6Avg(run,cnt) = meanMajOpininos6(end);
        cnt = cnt + 1;
    end
end

fig = figure;
scale = 2;
line_width = 2;
set(gcf, 'PaperPosition',[0.25 2.5 6 4]);
errorbar(list,mean(meanMajOpininos2Avg),std(meanMajOpininos2Avg)/sqrt(runs)*scale,'-k','linewidth',line_width);
hold on;
errorbar(list,mean(meanMajOpininos3Avg),std(meanMajOpininos3Avg)/sqrt(runs)*scale,'-.k','linewidth',line_width);
hold on;
errorbar(list,mean(meanMajOpininos5Avg),std(meanMajOpininos5Avg)/sqrt(runs)*scale,'-c','linewidth',line_width);
hold on;
errorbar(list,mean(meanMajOpininos6Avg),std(meanMajOpininos6Avg)/sqrt(runs)*scale,'--m','linewidth',line_width);
hold on;
errorbar(list,mean(meanMajOpininos4Avg),std(meanMajOpininos4Avg)/sqrt(runs)*scale,'-g','linewidth',line_width);
hold on;
errorbar(list,mean(meanMajOpininos1Avg),std(meanMajOpininos1Avg)/sqrt(runs)*scale,'-r','linewidth',line_width);
leg = legend('Degree','Betweenness','Closeness','Chen','EPN1','EPN2');
set(leg,'FontSize',7);
xlabel('\bf\alpha');
ylabel('\bfOpinion of population at The End of Simulation');
saveas(fig,'SocialPower.fig');

save('SocialData');


