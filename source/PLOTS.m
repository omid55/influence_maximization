close all;

fig = figure;
%plot(1:length(meanMajOpininos0Avg),meanMajOpininos0Avg,'-*',1:length(meanMajOpininos1Avg),meanMajOpininos1Avg,'-x',1:length(meanMajOpininos2Avg),meanMajOpininos2Avg,'-o',1:length(meanMajOpininos3Avg),meanMajOpininos3Avg,'-p');
plot(1:length(meanMajOpininos1Avg),meanMajOpininos1Avg,'-x',1:length(meanMajOpininos2Avg),meanMajOpininos2Avg,'-o',1:length(meanMajOpininos3Avg),meanMajOpininos3Avg,'-p',1:length(meanMajOpininos4Avg),meanMajOpininos4Avg,'-*',1:length(meanMajOpininos5Avg),meanMajOpininos5Avg,'-^');
xlabel('\bfTimestep');
ylabel('\bfAverage of Opinions');
set(gca,'XScale','log');
title('\bfTemporal Evolution of the Opinion Formation Process');
%legend('Natural Majority','Random Majority','Degree Targeted Majority','Proposed Targeted Majority','Location','Best');
legend('EPN2','Degree','Betweenness','EPN1','Closeness','Location','Best');
name = [path '/' f1 '/Majority'];
saveas(fig,[name '.fig']);
eval(['print -loos -dtiff ' name '.tiff']);

%% Last
fig = figure;
plot(sort(dist1), linspace(1-1/length(dist1), 0, length(dist1)),'-r','LineWidth',2);
hold on;
plot(sort(dist2), linspace(1-1/length(dist2), 0, length(dist2)),'-b','LineWidth',2);
hold on;
plot(sort(dist3), linspace(1-1/length(dist3), 0, length(dist3)),'-g','LineWidth',2);
hold on;
plot(sort(dist4), linspace(1-1/length(dist4), 0, length(dist4)),'-c','LineWidth',2);
hold on;
plot(sort(dist5), linspace(1-1/length(dist5), 0, length(dist5)),'-k','LineWidth',2);
legend('EPN2','Degree','Betweenness','Proposed1','Closeness');
ylabel('\bfP(X>o)');
xlabel('\bfOpinion');
title('\bfLast CCDF');
name = [path '/' f1 '/last_CCDF'];
saveas(fig,[name '.fig']);
eval(['print -loos -dtiff ' name '.tiff']);

fig = figure;
bin = -1:0.01:1;
x = (bin(1:length(bin)-1)+bin(2:length(bin)))/2;
plot(x,MyBinsCounting(dist1,bin),x,MyBinsCounting(dist2,bin),x,MyBinsCounting(dist3,bin),x,MyBinsCounting(dist4,bin),x,MyBinsCounting(dist5,bin),'LineWidth',2);
legend('EPN2','Degree','Betweenness','Proposed1','Closeness');
xlabel('\bfOpinion');
ylabel('\bfVolume');
title('\bfLast Dist');
name = [path '/' f1 '/last_Dist'];
saveas(fig,[name '.fig']);
eval(['print -loos -dtiff ' name '.tiff']);

%% Full
fig = figure;
plot(sort(fulldist1), linspace(1-1/length(fulldist1), 0, length(fulldist1)),'-r','LineWidth',2);
hold on;
plot(sort(fulldist2), linspace(1-1/length(fulldist2), 0, length(fulldist2)),'-b','LineWidth',2);
hold on;
plot(sort(fulldist3), linspace(1-1/length(fulldist3), 0, length(fulldist3)),'-g','LineWidth',2);
hold on;
plot(sort(fulldist4), linspace(1-1/length(fulldist4), 0, length(fulldist4)),'-c','LineWidth',2);
hold on;
plot(sort(fulldist5), linspace(1-1/length(fulldist5), 0, length(fulldist5)),'-k','LineWidth',2);
legend('EPN2','Degree','Betweenness','Proposed1','Closeness');
ylabel('\bfP(X>o)');
xlabel('\bfOpinion');
title('\bfFull CCDF');
name = [path '/' f1 '/full_CCDF'];
saveas(fig,[name '.fig']);
eval(['print -loos -dtiff ' name '.tiff']);

fig = figure;
bin = -1:0.01:1;
x = (bin(1:length(bin)-1)+bin(2:length(bin)))/2;
plot(x,MyBinsCounting(fulldist1,bin),x,MyBinsCounting(fulldist2,bin),x,MyBinsCounting(fulldist3,bin),x,MyBinsCounting(fulldist4,bin),x,MyBinsCounting(fulldist5,bin),'LineWidth',2);
legend('EPN2','Degree','Betweenness','Proposed1','Closeness');
xlabel('\bfOpinion');
ylabel('\bfVolume');
title('\bfFull Dist');
name = [path '/' f1 '/fullDist'];
saveas(fig,[name '.fig']);
eval(['print -loos -dtiff ' name '.tiff']);

fig = figure;
plot(1:length(meanFol1Avg),meanFol1Avg,'-x',1:length(meanFol2Avg),meanFol2Avg,'-o',1:length(meanFol3Avg),meanFol3Avg,'-p',1:length(meanFol4Avg),meanFol4Avg,'-*',1:length(meanFol5Avg),meanFol5Avg,'-^');
xlabel('\bfTimestep');
ylabel('\bfFollower Ratio');
set(gca,'XScale','log');
title('\bfTemporal Evolution of the Opinion Formation Process');
%legend('Natural Majority','Random Majority','Degree Targeted Majority','Proposed Targeted Majority','Location','Best');
legend('EPN2','Degree','Betweenness','EPN1','Closeness','Location','Best');
name = [path '/' f1 '/Followers'];
saveas(fig,[name '.fig']);
eval(['print -loos -dtiff ' name '.tiff']);