%Omid55
function [  ] = MuTunning(  )

%% Clear Everything
clc;
close all;


%% Initiallization
N = 1000;
InformedAgentsSize = ceil(N/10);
MaximumSimulationSteps = 10000;


%% Network Structure
averageDegree = 12;
%BA
A = BAnet(N,averageDegree/2,averageDegree/2);
sp = sparse(A);
net = CreateMap(sp);
f1 = 'BA';
disp(f1);


%% Simulation
runs = 10;
results = zeros(runs,10);
for run=1:runs
    fprintf(1,'Algorithm Run %d\n',run);
    
    %% Opinion initiallization
    % degree based
    scale = 0.7;
    degs = zeros(N,1);
    for i=1:N
        degs(i) = length(net(num2str(i)));
    end
    x = ((-1) .^ round(rand(N,1))) .* (scale * degs / max(degs));
        
    cnt = 1;
    for mu=0.1:0.1:1
        %fprintf(1,'Algorithm Run %d with mu=%d\n',run,mu);
        meanMajOpininos3 = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,3);
        results(run,cnt) = meanMajOpininos3(end);
        cnt = cnt + 1;
    end
end

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
errorbar((1:size(results,2))/10,mean(results),20*std(results)/sqrt(runs));
leg = legend('Last Opinion','Location','Northeast');
set(leg,'FontSize',7);
xlabel('\bf\mu');
ylabel('\bfOpinion of population at The End of Simulation');
%set(gca,'XScale','log');
%set(gca,'YScale','log');
saveas(fig, 'MuTunning.fig');
print -loos -dtiff MuTunning.tiff

save('MuTunningData.mat');


end

