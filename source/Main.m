%Omid55
%Complex Dynamical Networks' Mini Project
function [  ] = Main(  )

%% Clear Every thing
clc;
close all;
%clear all;


%% Initiallization
N = 100;
InformedAgentsSize = ceil(N/10);
MaximumSimulationSteps = 2000;
mu = 0.9;
muPrime = 0.9;


%% Simulation of the algorithm
XD = 1;
UI = 2;
alpha = 1;      % social power effection


%% Network Structure
averageDegree = 12;

%BA
A = BAnet(N,averageDegree/2,averageDegree/2);
sp = sparse(A);
net = CreateMap(sp);

% %ER
% erdosProbability = averageDegree/N;
% sp = erdos_reyni(N,erdosProbability);
% net = CreateMap(sp);
% %erdosProbability = averageDegree/2*N;
% %[net,sp] = ErdosRenyi(N,erdosProbability);

% %Regular
% sp = CreateRegularLattice(N,averageDegree);
% net = CreateMap(sp);

% % %WS
% % sp = WattsStrogatzCreator(N,averageDegree,0.1);
% % net = CreateMap(sp);
% % REAL NETWORK
% network_file = 'out.ego-facebook-friendships-net';
% a = importdata(network_file);
% MM = max(max(a(:,1)),max(a(:,2)));
% sp = sparse(a(:,1),a(:,2),1,MM,MM);
% net = CreateMap(sp);
% N = size(sp,1);
% InformedAgentsSize = ceil(N/10);
% %MaximumSimulationSteps = 5000;
% f1 = 'FacebookFriendships';

% %%%
% %my own network
% ed = [1 2;1 3;1 4;1 5;1 7;5 6;4 6;6 7;7 8;9 8;9 13;9 10;9 11;9 12;11 12;10 11;10 12];
% N = max(max(ed));
% sp = sparse(ed(:,1),ed(:,2),1,N,N);
% sp = max(sp,sp');
% net = CreateMap(sp);
% InformedAgentsSize = 3;
% %%%


%% Simulation
runs = 2;
leng = MaximumSimulationSteps + 1;
meanAgOpininos0Avg = zeros(1,leng);
meanMajOpininos0Avg = zeros(1,leng);
meanAgOpininos1Avg = zeros(1,leng);
meanMajOpininos1Avg = zeros(1,leng);
meanAgOpininos2Avg = zeros(1,leng);
meanMajOpininos2Avg = zeros(1,leng);
meanAgOpininos3Avg = zeros(1,leng);
meanMajOpininos3Avg = zeros(1,leng);

meanAllOpininos0AvgMat = zeros(runs,leng);
meanAllOpininos1AvgMat = zeros(runs,leng);
meanAllOpininos2AvgMat = zeros(runs,leng);
meanAllOpininos3AvgMat = zeros(runs,leng);

for run=1:runs

    fprintf(1,'Algorithm Run %d\n',run);
    
    net0 = containers.Map(net.keys,net.values);
    net1 = containers.Map(net.keys,net.values);
    net2 = containers.Map(net.keys,net.values);
    net3 = containers.Map(net.keys,net.values);

   %% Opinion initiallization
%     % by random
%     x = 2 * rand(N+InformedAgentsSize,1) - 1;

    % degree based
    scale = 1;
    %degs = sum(sp);
    degs = zeros(N,1);
    for i=1:N
        degs(i) = length(net(num2str(i)));
    end
    %x = ((-1) .^ round(rand(N,1))) .* (scale * degs / max(degs));
    x = ((-1) .^ round(rand(N,1))) .* power(degs/max(degs),1/2);
    
    [meanAllOpinions0,meanAgOpininos0,meanMajOpininos0] = SimulationMethod(x(1:N),N,0,MaximumSimulationSteps,mu,XD,UI,alpha,sp,net0,1);
    [meanAllOpinions1,meanAgOpininos1,meanMajOpininos1] = SimulationMethod(x,N,InformedAgentsSize,MaximumSimulationSteps,mu,XD,UI,alpha,sp,net1,1);
    [meanAllOpinions2,meanAgOpininos2,meanMajOpininos2] = SimulationMethod(x,N,InformedAgentsSize,MaximumSimulationSteps,mu,XD,UI,alpha,sp,net2,2);
    [meanAllOpinions3,meanAgOpininos3,meanMajOpininos3] = SimulationMethod(x,N,InformedAgentsSize,MaximumSimulationSteps,mu,XD,UI,alpha,sp,net3,3);
    
    meanAllOpininos0AvgMat(run,:) = meanAllOpinions0;
    meanAllOpininos1AvgMat(run,:) = meanAllOpinions1;
    meanAllOpininos2AvgMat(run,:) = meanAllOpinions2;
    meanAllOpininos3AvgMat(run,:) = meanAllOpinions3;
    
    meanAgOpininos0Avg = meanAgOpininos0Avg + meanAgOpininos0;
    meanMajOpininos0Avg = meanMajOpininos0Avg + meanMajOpininos0;
    meanAgOpininos1Avg = meanAgOpininos1Avg + meanAgOpininos1;
    meanMajOpininos1Avg = meanMajOpininos1Avg + meanMajOpininos1;
    meanAgOpininos2Avg = meanAgOpininos2Avg + meanAgOpininos2;
    meanMajOpininos2Avg = meanMajOpininos2Avg + meanMajOpininos2;
    meanAgOpininos3Avg = meanAgOpininos3Avg + meanAgOpininos3;
    meanMajOpininos3Avg = meanMajOpininos3Avg + meanMajOpininos3;
    
end

meanAgOpininos0Avg = meanAgOpininos0Avg / runs;
meanMajOpininos0Avg = meanMajOpininos0Avg / runs;
meanAgOpininos1Avg = meanAgOpininos1Avg / runs;
meanMajOpininos1Avg = meanMajOpininos1Avg / runs;
meanAgOpininos2Avg = meanAgOpininos2Avg / runs;
meanMajOpininos2Avg = meanMajOpininos2Avg / runs;
meanAgOpininos3Avg = meanAgOpininos3Avg / runs;
meanMajOpininos3Avg = meanMajOpininos3Avg / runs;
    

fig = figure;
plot(1:leng,mean(meanAllOpininos0AvgMat),1:leng,mean(meanAllOpininos1AvgMat),1:leng,mean(meanAllOpininos2AvgMat),1:leng,mean(meanAllOpininos3AvgMat));
xlabel('Timestep');
ylabel('Average of Opinions');
title('All Society Temporal Evolution of the Opinion Formation Process');
legend('Natural Whole Society','Random Whole Society','Degree Targeted Whole Society','Proposed Targeted Whole Society','Location','Best');
saveas(fig,'All.fig');
print -loos -dtiff All.tiff;
 
fig = figure;
plot(1:length(meanAgOpininos0Avg),meanAgOpininos0Avg,1:length(meanAgOpininos1Avg),meanAgOpininos1Avg,1:length(meanMajOpininos1Avg),meanMajOpininos1Avg,1:length(meanAgOpininos2Avg),meanAgOpininos2Avg,1:length(meanMajOpininos2Avg),meanMajOpininos2Avg,1:length(meanAgOpininos3Avg),meanAgOpininos3Avg,1:length(meanMajOpininos3Avg),meanMajOpininos3Avg);
xlabel('Timestep');
ylabel('Average of Opinions');
title('Temporal Evolution of the Opinion Formation Process');
legend('Natural Informed Agents','Random Informed Agents','Random Majority','Degree Targeted Informed Agents','Degree Targeted Majority','Proposed Targeted Informed Agents','Proposed Targeted Majority','Location','Best');
saveas(fig,'Separated.fig');
print -loos -dtiff Separated.tiff;

fig = figure;
plot(1:length(meanMajOpininos0Avg),meanMajOpininos0Avg,1:length(meanMajOpininos1Avg),meanMajOpininos1Avg,1:length(meanMajOpininos2Avg),meanMajOpininos2Avg,1:length(meanMajOpininos3Avg),meanMajOpininos3Avg);
xlabel('Timestep');
ylabel('Average of Opinions');
title('Temporal Evolution of the Opinion Formation Process');
legend('Natural Majority','Random Majority','Degree Targeted Majority','Proposed Targeted Majority','Location','Best');
saveas(fig,'Majority.fig');
print -loos -dtiff Majority.tiff;

fig = figure;
hold on;
set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
errorbar(1:size(meanAllOpininos1AvgMat,2),mean(meanAllOpininos1AvgMat),std(meanAllOpininos1AvgMat),'r')
errorbar(1:size(meanAllOpininos2AvgMat,2),mean(meanAllOpininos2AvgMat),std(meanAllOpininos2AvgMat),'b')
errorbar(1:size(meanAllOpininos3AvgMat,2),mean(meanAllOpininos3AvgMat),std(meanAllOpininos3AvgMat),'g')
leg = legend('Random Majority','Degree Targeted Majority','Proposed Targeted Majority','Location','Northeast');
set(leg,'FontSize',7);
xlabel('Timestep');
ylabel('Average of Opinions');
set(gca,'XScale','log');
%set(gca,'YScale','log');
saveas(fig,'ErrorBar.fig');
print -loos -dtiff ErrorBar.tiff;



save('FullData.mat');

end

