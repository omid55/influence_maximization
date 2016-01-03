%Omid55
function [  ] = RunMain( path,netType )

%% Initiallization
N = 100;
InformedAgentsSize = ceil(N/10);
MaximumSimulationSteps = 1000;
mu = 0.8;
runs = 2;
Randomness = 0.01;
scale = 0.7;


%% Network Structure
averageDegree = 12;

switch netType
    
    case 1
    %BA
    A = BAnet(N,averageDegree/2,averageDegree/2);
    sp = sparse(A);
    %net = CreateMap(sp);
    f1 = 'BA';
    disp(f1);

    case 2
    %ER
    erdosProbability = averageDegree/N;
    sp = erdos_reyni(N,erdosProbability);
    %net = CreateMap(sp);
    %erdosProbability = averageDegree/2*N;
    %[net,sp] = ErdosRenyi(N,erdosProbability);
    f1 = 'ER';
    disp(f1);

    case 3
    %Modular BA
    coef = 5;
    probability = 0.1;
    A = ModularBA(coef,N/coef,averageDegree/2,averageDegree/2,probability);
    sp = sparse(A);
    %net = CreateMap(sp);
    f1 = 'ModularBA';
    disp(f1);

    case 4
    %WS
    sp = WattsStrogatzCreator(N,averageDegree,0.1);
    %net = CreateMap(sp);
    f1 = 'WS';
    disp(f1);

    case 5
    %Modular WS
    coef = 5;
    probability = 0.1;
    p_WS = 0.1;
    A = ModularWS(coef,N/coef,averageDegree,p_WS,probability);
    sp = sparse(A);
    %net = CreateMap(sp);
    f1 = 'ModularWS';
    disp(f1);
    
    case 6
    %Forest Fire
    sp = ForestFireCreator(N,0.2,0.2);
    f1 = 'FF';
    disp(f1);
    
    case 7
    % REAL NETWORK
    network_file = 'out - ucidata-zachary.txt';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM);
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'Zachary';
    %runs = 10;
    disp(f1);    
    
    case 8
    network_file = 'out - advogato.txt';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM);
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 2 * N;
    f1 = 'Advogato';
    %runs = 10;
    disp(f1);
    
    case 9
    network_file = 'out - elec.txt';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM);
    sp = max(sp,sp');
    sp(find(sp>1)) = 1;
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'Elec';
    %runs = 10;
    disp(f1);
    
    case 10
    network_file = 'out - petster-hamster.txt';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM);
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 10000;
    f1 = 'Hamster';
    %runs = 10;
    disp(f1);
    
    case 11
    network_file = 'out.opsahl-ucforum';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM);
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'UCIrvine';
    %runs = 10;
    disp(f1);
    
    case 12
    network_file = 'out.petster-friendships-hamster';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM);
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'HamsterFriendships';
    %runs = 10;
    disp(f1);
    
% %     case 11
% %     network_file = 'out - ego-facebook.txt';
% %     a = importdata(network_file);
% %     MM = max(max(a(:,1)),max(a(:,2)));
% %     sp = sparse(a(:,1),a(:,2),1,MM,MM);
% %     sp = max(sp,sp');
% %     %net = CreateMap(sp);
% %     N = size(sp,1);
% %     InformedAgentsSize = ceil(N/10);
% %     %MaximumSimulationSteps = 5000;
% %     f1 = 'Facebook';
% %     %runs = 10;
% %     disp(f1);
% %     
% %     case 12
% %     network_file = 'facebook-wall.txt.anon';
% %     a = importdata(network_file);
% %     MM = max(max(a(:,1)),max(a(:,2)));
% %     sp = sparse(a(:,1),a(:,2),1,MM,MM);
% %     sp = max(sp,sp');
% %     %net = CreateMap(sp);
% %     N = size(sp,1);
% %     InformedAgentsSize = ceil(N/10);
% %     %MaximumSimulationSteps = 5000;
% %     f1 = 'FacebookWallLarge';
% %     %runs = 10;
% %     disp(f1);
    
    case 13
    %%%
    %my own network
    ed = [1 2;1 3;1 4;1 5;1 7;5 6;4 6;6 7;7 8;9 8;9 13;9 10;9 11;9 12;11 12;10 11;10 12];
    N = max(max(ed));
    sp = sparse(ed(:,1),ed(:,2),1,N,N);
    sp = max(sp,sp');
    %net = CreateMap(sp);
    InformedAgentsSize = ceil(N/10)+1;
    f1 = 'Sample';
    disp(f1);
    %%%

end


%% Simulation
leng = MaximumSimulationSteps + 1;
meanMajOpininos1Avg = zeros(leng,1);
meanMajOpininos2Avg = zeros(leng,1);
meanMajOpininos3Avg = zeros(leng,1);
meanMajOpininos4Avg = zeros(leng,1);

% %Adjacent Calculation
A = full(sp);
degs = sum(A)';

for run=1:runs

    fprintf(1,'Algorithm Run %d\n',run);
    
   %% Opinion initiallization
    % degree based
    x = ((-1) .^ round(rand(N,1))) .* (scale * degs / max(degs));
    %%x = ((-1) .^ round(rand(N,1))) .* power(degs/max(degs),1/2);
    
%     % random based
%     x = -1 + 2 * rand(N,1);
    
    % Adjacent Selection Approaches
    f2 = 'Single';

    meanMajOpininos1 = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,1,Randomness)';
    meanMajOpininos2 = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,2,Randomness)';
    meanMajOpininos3 = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,3,Randomness)';
    meanMajOpininos4 = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,4,Randomness)';
    
    %meanMajOpininos0Avg = meanMajOpininos0Avg + meanMajOpininos0';
    meanMajOpininos1Avg = meanMajOpininos1Avg + meanMajOpininos1';
    meanMajOpininos2Avg = meanMajOpininos2Avg + meanMajOpininos2';
    meanMajOpininos3Avg = meanMajOpininos3Avg + meanMajOpininos3';
    meanMajOpininos4Avg = meanMajOpininos4Avg + meanMajOpininos4';
end

%meanMajOpininos0Avg = meanMajOpininos0Avg / runs;
meanMajOpininos1Avg = meanMajOpininos1Avg / runs;
meanMajOpininos2Avg = meanMajOpininos2Avg / runs;
meanMajOpininos3Avg = meanMajOpininos3Avg / runs;
meanMajOpininos4Avg = meanMajOpininos4Avg / runs;
 

%% Results
mkdir([path '/' f1 '/' f2]);

fig1 = figure;
%plot(1:length(meanMajOpininos0Avg),meanMajOpininos0Avg,'-*',1:length(meanMajOpininos1Avg),meanMajOpininos1Avg,'-x',1:length(meanMajOpininos2Avg),meanMajOpininos2Avg,'-o',1:length(meanMajOpininos3Avg),meanMajOpininos3Avg,'-p');
plot(1:length(meanMajOpininos1Avg),meanMajOpininos1Avg,'-x',1:length(meanMajOpininos2Avg),meanMajOpininos2Avg,'-o',1:length(meanMajOpininos3Avg),meanMajOpininos3Avg,'-p',1:length(meanMajOpininos4Avg),meanMajOpininos4Avg,'-*');
xlabel('Timestep');
ylabel('Average of Opinions');
%set(gca,'XScale','log');
title('Temporal Evolution of the Opinion Formation Process');
%legend('Natural Majority','Random Majority','Degree Targeted Majority','Proposed Targeted Majority','Location','Best');
legend('Random','Degree','Inverse Degree','Proposed','Location','Best');
name = [path '/' f1 '/' f2 '/Majority'];
saveas(fig1,[name '.fig']);
eval(['print -loos -dtiff ' name '.tiff']);
eval(['print -loos -dtiff ' [f1 '_' f2 '_Majority'] '.tiff']); 

save([path '/' f1 '/' f2 '/FullData.mat']);

close all;
clc;
pause(3);  % for closing all of figures

end

