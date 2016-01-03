%% PLOTTING
close all;
COUNT = 20;
set(gcf, 'PaperPosition',[0.25 2.5 12 6]);
figure(1);
width = 2;
len = length(meanMajOpininos2Avg);

hold on;
indices = 1:floor(len/COUNT):len;
m2 = mean(meanMajOpininos2Avg);
m3 = mean(meanMajOpininos3Avg);
m5 = mean(meanMajOpininos5Avg);
m4 = mean(meanMajOpininos4Avg);
m1 = mean(meanMajOpininos1Avg);

s2 = std(meanMajOpininos2Avg)/sqrt(runs);
s3 = std(meanMajOpininos3Avg)/sqrt(runs);
s5 = std(meanMajOpininos5Avg)/sqrt(runs);
s4 = std(meanMajOpininos4Avg)/sqrt(runs);
s1 = std(meanMajOpininos1Avg)/sqrt(runs);

errorbar(m2(indices),s2(indices),'-k','linewidth',width);
errorbar(m3(indices),s3(indices),'-.k','linewidth',width);
errorbar(m5(indices),s5(indices),'-c','linewidth',width);
errorbar(m4(indices),s4(indices),'--c','linewidth',width);
errorbar(m1(indices),s1(indices),'-r','linewidth',width);

axis([1 12 0 1]);
xlabel('\bfTimestep','FontSize',12);
legend('Degree','Betweenness','Closeness','EPN1','EPN2','Location','SouthEast');
ylabel('\bfAverage of Opinions','FontSize',12);
%% PLOTTING
