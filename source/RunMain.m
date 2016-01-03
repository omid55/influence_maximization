%Omid55
function [  ] = RunMain( path,netType,type )

%% Initiallization
% addpath('/agbs/cluster/oaskaris/im/mb/matlab_bgl');
NET_TYPE = netType
N = 1000;
InformedAgentsSize = ceil(N/10);

if type == 1
    % FIG5A
    MaximumSimulationSteps = 10000;
    mu = 0.4;
    alpha = 0;
else
    if type == 2
        % FIG5B
        MaximumSimulationSteps = 20000;
        mu = 0.8;
        alpha = 1;
    else
        % REAL
        MaximumSimulationSteps = 40000;
        mu = 0.2;
        alpha = 0;
    end
end

runs = 10;
Randomness = 0.01;
beta = 0.1;


% NET_TYPE = netType
% N = 2000;
% InformedAgentsSize = ceil(N/10);
% MaximumSimulationSteps = 5000;
% mu = 0.7; 
% runs = 100;
% Randomness = 0.01;
% alpha = 0;
% beta = 0.1;


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
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'Zachary';
    disp(f1);    
    
    case 8
    network_file = 'out - advogato.txt';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 2 * N;
    f1 = 'Advogato';
    disp(f1);
    
    case 9
    network_file = 'out - elec.txt';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    sp(find(sp>1)) = 1;
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'Elec';
    disp(f1);
    
    case 10
    network_file = 'out - petster-hamster.txt';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 10000;
    f1 = 'Hamster';
    disp(f1);
    
    case 11
    network_file = 'out.petster-friendships-hamster';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
%     sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'HamsterFriendships';
    disp(f1);
	
    case 12
    network_file = 'out - ego-facebook.txt';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'Facebook';
    disp(f1);

    case 13
    network_file = 'opsahl';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
%     sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'UCIrvine';
    disp(f1);
    
    case 14
    network_file = 'out.facebook-wosn-links';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'Facebookfriendships';
    disp(f1);
       
    case 15
    network_file = 'out.ego-gplus.txt';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'GooglePlus';
    disp(f1);
	
    case 17
    network_file = 'out.loc-brightkite_edges.txt';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'Brightkite';
    disp(f1);
    
    case 18
    network_file = 'out.web-Google';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'WebGoogle';
    disp(f1);
    
    case 16
    network_file = 'out.petster-friendships-dog';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'PetsterFriendships';
    disp(f1);
    
    case 19
    network_file = 'out.orkut-links';
    a = importdata(network_file);
    MM = max(max(a(:,1)),max(a(:,2)));
    sp = sparse(a(:,1),a(:,2),1,MM,MM); clear('a');
    sp = max(sp,sp');
    %net = CreateMap(sp);
    N = size(sp,1);
    InformedAgentsSize = ceil(N/10);
    %MaximumSimulationSteps = 5000;
    f1 = 'Orkut';
    disp(f1);
    
    
%     case 13
%     %%%
%     %my own network
%     ed = [1 2;1 3;1 4;1 5;1 7;5 6;4 6;6 7;7 8;9 8;9 13;9 10;9 11;9 12;11 12;10 11;10 12];
%     N = max(max(ed));
%     sp = sparse(ed(:,1),ed(:,2),1,N,N);
%     sp = max(sp,sp');
%     %net = CreateMap(sp);
%     InformedAgentsSize = ceil(N/10)+1;
%     f1 = 'Sample';
%     disp(f1);
%     %%%
    
end


%% Simulation
leng = MaximumSimulationSteps + 1;
meanMajOpininos1Avg = zeros(runs,leng);
% meanMajOpininos2Avg = zeros(leng,1);
% meanMajOpininos3Avg = zeros(leng,1);
% meanMajOpininos4Avg = zeros(leng,1);
% meanMajOpininos5Avg = zeros(leng,1);

meanMajOpininos1Avg = zeros(runs,leng);
meanMajOpininos2Avg = zeros(runs,leng);
meanMajOpininos3Avg = zeros(runs,leng);
meanMajOpininos4Avg = zeros(runs,leng);
meanMajOpininos5Avg = zeros(runs,leng);
meanMajOpininos6Avg = zeros(runs,leng);

meanFol1Avg = zeros(leng,1);
meanFol2Avg = zeros(leng,1);
meanFol3Avg = zeros(leng,1);
meanFol4Avg = zeros(leng,1);
meanFol5Avg = zeros(leng,1);
meanFol6Avg = zeros(leng,1);

%% Betweenness & Closeness Calculation
betweennesses = betweenness_centrality(sp);
[D ~] = all_shortest_paths(sp,struct('algname','floyd_warshall'));
closenesses = (N-1) ./ sum(D);

% %Adjacent Calculation
A = full(sp);
clear('sp');
% degs = sum(A)';
dist1 = zeros(1,N+InformedAgentsSize);
dist2 = zeros(1,N+InformedAgentsSize);
dist3 = zeros(1,N+InformedAgentsSize);
dist4 = zeros(1,N+InformedAgentsSize);
dist5 = zeros(1,N+InformedAgentsSize);
dist6 = zeros(1,N+InformedAgentsSize);
fulldist1 = zeros(1,N+InformedAgentsSize);
fulldist2 = zeros(1,N+InformedAgentsSize);
fulldist3 = zeros(1,N+InformedAgentsSize);
fulldist4 = zeros(1,N+InformedAgentsSize);
fulldist5 = zeros(1,N+InformedAgentsSize);
fulldist6 = zeros(1,N+InformedAgentsSize);

run = 1;
while run <= runs

    fprintf(1,'Algorithm Run %d\n',run);
    x = -1 + 2 * rand(N,1);

    [meanMajOpininos1,d1,fd1,fol1] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,1,Randomness,alpha,beta,betweennesses,closenesses);    
    [meanMajOpininos2,d2,fd2,fol2] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,2,Randomness,alpha,beta,betweennesses,closenesses);
    [meanMajOpininos3,d3,fd3,fol3] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,3,Randomness,alpha,beta,betweennesses,closenesses);
    [meanMajOpininos4,d4,fd4,fol4] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,4,Randomness,alpha,beta,betweennesses,closenesses);
    [meanMajOpininos5,d5,fd5,fol5] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,5,Randomness,alpha,beta,betweennesses,closenesses);
    [meanMajOpininos6,d6,fd6,fol6] = SimulationMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,6,Randomness,alpha,beta,betweennesses,closenesses);
    
    dist1 = dist1 + d1/runs;
    dist2 = dist2 + d2/runs;
    dist3 = dist3 + d3/runs;
    dist4 = dist4 + d4/runs;
    dist5 = dist5 + d5/runs;
    dist6 = dist6 + d6/runs;
    
    fulldist1 = fulldist1 + fd1/runs;
    fulldist2 = fulldist2 + fd2/runs;
    fulldist3 = fulldist3 + fd3/runs;
    fulldist4 = fulldist4 + fd4/runs;
    fulldist5 = fulldist5 + fd5/runs;
    fulldist6 = fulldist6 + fd6/runs;
    
%     meanMajOpininos1Avg = meanMajOpininos1Avg + meanMajOpininos1'/runs;
%     meanMajOpininos2Avg = meanMajOpininos2Avg + meanMajOpininos2'/runs;
%     meanMajOpininos3Avg = meanMajOpininos3Avg + meanMajOpininos3'/runs;
%     meanMajOpininos4Avg = meanMajOpininos4Avg + meanMajOpininos4'/runs;
%     meanMajOpininos5Avg = meanMajOpininos5Avg + meanMajOpininos5'/runs;
%     
    meanMajOpininos1Avg(run,:) = meanMajOpininos1;
    meanMajOpininos2Avg(run,:) = meanMajOpininos2;
    meanMajOpininos3Avg(run,:) = meanMajOpininos3;
    meanMajOpininos4Avg(run,:) = meanMajOpininos4;
    meanMajOpininos5Avg(run,:) = meanMajOpininos5;
    meanMajOpininos6Avg(run,:) = meanMajOpininos6;
    
    meanFol1Avg = meanFol1Avg + fol1'/runs;
    meanFol2Avg = meanFol2Avg + fol2'/runs;
    meanFol3Avg = meanFol3Avg + fol3'/runs;
    meanFol4Avg = meanFol4Avg + fol4'/runs;
    meanFol5Avg = meanFol5Avg + fol5'/runs;
    meanFol6Avg = meanFol6Avg + fol6'/runs;
    
    run = run + 1;
    
end


%% Results
mkdir([path '/' f1]);
save([path '/' f1 '/FullData.mat']);


%% PLOTTING
COUNT = 20;    %how many points you need in last figure
close all;
set(gcf, 'PaperPosition',[0.25 2.5 12 6]);
fig1 = figure(1);
width = 2;
len = length(meanMajOpininos2Avg);

hold on;
indices = 1:floor(len/COUNT):len;
m2 = mean(meanMajOpininos2Avg);
m3 = mean(meanMajOpininos3Avg);
m5 = mean(meanMajOpininos5Avg);
m4 = mean(meanMajOpininos4Avg);
m1 = mean(meanMajOpininos1Avg);
m6 = mean(meanMajOpininos6Avg);

scale = 4;
s2 = scale * std(meanMajOpininos2Avg);    %std(meanMajOpininos2Avg)/sqrt(runs);
s3 = scale * std(meanMajOpininos3Avg);    %std(meanMajOpininos3Avg)/sqrt(runs);
s5 = scale * std(meanMajOpininos5Avg);    %std(meanMajOpininos5Avg)/sqrt(runs);
s4 = scale * std(meanMajOpininos4Avg);    %std(meanMajOpininos4Avg)/sqrt(runs);
s1 = scale * std(meanMajOpininos1Avg);     %std(meanMajOpininos1Avg)/sqrt(runs);
s6 = scale * std(meanMajOpininos6Avg);    %std(meanMajOpininos6Avg)/sqrt(runs);


errorbar(m2(indices),s2(indices),'-k','linewidth',width);
errorbar(m3(indices),s3(indices),'-.k','linewidth',width);
errorbar(m5(indices),s5(indices),'-c','linewidth',width);
errorbar(m6(indices),s6(indices),'--m','linewidth',width);
errorbar(m4(indices),s4(indices),'-g','linewidth',width);
errorbar(m1(indices),s1(indices),'-r','linewidth',width);

xlabel('\bfTimestep','FontSize',12);
legend('Degree','Betweenness','Closeness','Chen','EPN1','EPN2','Location','SouthEast');
ylabel('\bfAverage of Opinions','FontSize',12);

axis([1 15 0 1]);
set(gca,'XTickLabel',{'3000','6000','9000','12000','15000'},'XTick',[3 6 9 12 15]);
%% PLOTTING


% close all;

% fig = figure;
% %plot(1:length(meanMajOpininos0Avg),meanMajOpininos0Avg,'-*',1:length(meanMajOpininos1Avg),meanMajOpininos1Avg,'-x',1:length(meanMajOpininos2Avg),meanMajOpininos2Avg,'-o',1:length(meanMajOpininos3Avg),meanMajOpininos3Avg,'-p');
% plot(1:length(meanMajOpininos1Avg),meanMajOpininos1Avg,'-x',1:length(meanMajOpininos2Avg),meanMajOpininos2Avg,'-o',1:length(meanMajOpininos3Avg),meanMajOpininos3Avg,'-p',1:length(meanMajOpininos4Avg),meanMajOpininos4Avg,'-*',1:length(meanMajOpininos5Avg),meanMajOpininos5Avg,'-^');
% xlabel('\bfTimestep');
% ylabel('\bfAverage of Opinions');
% set(gca,'XScale','log');
% title('\bfTemporal Evolution of the Opinion Formation Process');
% %legend('Natural Majority','Random Majority','Degree Targeted Majority','Proposed Targeted Majority','Location','Best');
% legend('EPN2','Degree','Betweenness','EPN1','Closeness','Location','Best');
% name = [path '/' f1 '/Majority'];
% saveas(fig,[name '.fig']);
% eval(['print -loos -dtiff ' name '.tiff']);

% %% Last
% fig = figure;
% plot(sort(dist1), linspace(1-1/length(dist1), 0, length(dist1)),'-r','LineWidth',2);
% hold on;
% plot(sort(dist2), linspace(1-1/length(dist2), 0, length(dist2)),'-b','LineWidth',2);
% hold on;
% plot(sort(dist3), linspace(1-1/length(dist3), 0, length(dist3)),'-g','LineWidth',2);
% hold on;
% plot(sort(dist4), linspace(1-1/length(dist4), 0, length(dist4)),'-c','LineWidth',2);
% hold on;
% plot(sort(dist5), linspace(1-1/length(dist5), 0, length(dist5)),'-k','LineWidth',2);
% legend('EPN2','Degree','Betweenness','EPN1','Closeness');
% ylabel('\bfP(X>o)');
% xlabel('\bfOpinion');
% title('\bfLast CCDF');
% name = [path '/' f1 '/last_CCDF'];
% saveas(fig,[name '.fig']);
% eval(['print -loos -dtiff ' name '.tiff']);
% 
% fig = figure;
% bin = -1:0.01:1;
% x = (bin(1:length(bin)-1)+bin(2:length(bin)))/2;
% plot(x,MyBinsCounting(dist1,bin),x,MyBinsCounting(dist2,bin),x,MyBinsCounting(dist3,bin),x,MyBinsCounting(dist4,bin),x,MyBinsCounting(dist5,bin),'LineWidth',2);
% legend('EPN2','Degree','Betweenness','EPN1','Closeness');
% xlabel('\bfOpinion');
% ylabel('\bfVolume');
% title('\bfLast Dist');
% name = [path '/' f1 '/last_Dist'];
% saveas(fig,[name '.fig']);
% eval(['print -loos -dtiff ' name '.tiff']);
% 
% %% Full
% fig = figure;
% plot(sort(fulldist1), linspace(1-1/length(fulldist1), 0, length(fulldist1)),'-r','LineWidth',2);
% hold on;
% plot(sort(fulldist2), linspace(1-1/length(fulldist2), 0, length(fulldist2)),'-b','LineWidth',2);
% hold on;
% plot(sort(fulldist3), linspace(1-1/length(fulldist3), 0, length(fulldist3)),'-g','LineWidth',2);
% hold on;
% plot(sort(fulldist4), linspace(1-1/length(fulldist4), 0, length(fulldist4)),'-c','LineWidth',2);
% hold on;
% plot(sort(fulldist5), linspace(1-1/length(fulldist5), 0, length(fulldist5)),'-k','LineWidth',2);
% legend('EPN2','Degree','Betweenness','EPN1','Closeness');
% ylabel('\bfP(X>o)');
% xlabel('\bfOpinion');
% title('\bfFull CCDF');
% name = [path '/' f1 '/full_CCDF'];
% saveas(fig,[name '.fig']);
% eval(['print -loos -dtiff ' name '.tiff']);
% 
% fig = figure;
% bin = -1:0.01:1;
% x = (bin(1:length(bin)-1)+bin(2:length(bin)))/2;
% plot(x,MyBinsCounting(fulldist1,bin),x,MyBinsCounting(fulldist2,bin),x,MyBinsCounting(fulldist3,bin),x,MyBinsCounting(fulldist4,bin),x,MyBinsCounting(fulldist5,bin),'LineWidth',2);
% legend('EPN2','Degree','Betweenness','EPN1','Closeness');
% xlabel('\bfOpinion');
% ylabel('\bfVolume');
% title('\bfFull Dist');
% name = [path '/' f1 '/fullDist'];
% saveas(fig,[name '.fig']);
% eval(['print -loos -dtiff ' name '.tiff']);
% 
% fig = figure;
% plot(1:length(meanFol1Avg),meanFol1Avg,'-x',1:length(meanFol2Avg),meanFol2Avg,'-o',1:length(meanFol3Avg),meanFol3Avg,'-p',1:length(meanFol4Avg),meanFol4Avg,'-*',1:length(meanFol5Avg),meanFol5Avg,'-^');
% xlabel('\bfTimestep');
% ylabel('\bfFollower Ratio');
% set(gca,'XScale','log');
% title('\bfTemporal Evolution of the Opinion Formation Process');
% %legend('Natural Majority','Random Majority','Degree Targeted Majority','Proposed Targeted Majority','Location','Best');
% legend('EPN2','Degree','Betweenness','EPN1','Closeness','Location','Best');
% name = [path '/' f1 '/Followers'];
% saveas(fig,[name '.fig']);
% eval(['print -loos -dtiff ' name '.tiff']);


end

