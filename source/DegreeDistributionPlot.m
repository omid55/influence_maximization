%% Different Degree Distribution Ploting
close all;
clc;
clear;

N = 100;
InformedAgentsSize = ceil(N/10);
MaximumSimulationSteps = 10000;
mu = 0.9;
runs = 30;
Randomness = 0.01;
scale = 0.7;
alpha = 0;
path = 'Results';


%% Simulation
leng = MaximumSimulationSteps + 1;
val2Lev = zeros(10,1);
val1Lev = zeros(10,1);

for run=1:runs
    fprintf(1,'Algorithm Run %d\n',run);
    for averageDegree = 8 : 2 : 26
        A = BAnet(N,averageDegree/2,averageDegree/2);
        sp = sparse(A);

        x = -1 + 2 * rand(N,1);
        [meanMajOpininos1Lev,~,~,fol1] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,1,Randomness,alpha);
        [meanMajOpininos4Lev,~,~,fol4] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,4,Randomness,alpha);
        
        val2Lev((averageDegree-6)/2) = val2Lev((averageDegree-6)/2) + fol1(end) / runs;
        val1Lev((averageDegree-6)/2) = val1Lev((averageDegree-6)/2) + fol4(end) / runs;
    end
end

set(gcf, 'PaperPosition',[0.25 2.5 6 4]);
plot(8 : 2 : 26,val1Lev,'-x', 8 : 2 : 26,val2Lev,'-o');
legend('1 Level','2 Level');
print -dtiff -loose AverageDegree.tiff



