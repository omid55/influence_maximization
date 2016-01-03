% Omid55
% Result Generator Item

if ttype == 1
    
    %% Artificial Networks
    disp('Artificial Networks');
    set(gcf, 'PaperPosition',[0.25 2.5 12 6]);
    figure(1);
    subplot(2,3,my_cnt);
    % runs = 1;
    line_width = 2;
    len = length(meanMajOpininos2Avg);
    % plot(1:length(meanMajOpininos2Avg),meanMajOpininos2Avg,'-k',1:length(meanMajOpininos3Avg),meanMajOpininos3Avg,'-.k',1:length(meanMajOpininos5Avg),meanMajOpininos5Avg,'-c',1:length(meanMajOpininos4Avg),meanMajOpininos4Avg,'--c',1:length(meanMajOpininos1Avg),meanMajOpininos1Avg,'-r','LineWidth',3);

    hold on;
    indices = 1:floor(len/COUNT):len;
    m2 = mean(meanMajOpininos2Avg);
    m3 = mean(meanMajOpininos3Avg);
    m5 = mean(meanMajOpininos5Avg);
    m4 = mean(meanMajOpininos4Avg);
    m1 = mean(meanMajOpininos1Avg);
    m6 = mean(meanMajOpininos6Avg);

    scale = 3;
    s2 = scale * std(meanMajOpininos2Avg);    %std(meanMajOpininos2Avg)/sqrt(runs);
    s3 = scale * std(meanMajOpininos3Avg);    %std(meanMajOpininos3Avg)/sqrt(runs);
    s5 = scale * std(meanMajOpininos5Avg);    %std(meanMajOpininos5Avg)/sqrt(runs);
    s4 = scale * std(meanMajOpininos4Avg);    %std(meanMajOpininos4Avg)/sqrt(runs);
    s1 = scale * std(meanMajOpininos1Avg);     %std(meanMajOpininos1Avg)/sqrt(runs);
    s6 = scale * std(meanMajOpininos6Avg);    %std(meanMajOpininos6Avg)/sqrt(runs);

    errorbar(m2(indices),s2(indices),'-k','linewidth',line_width);
    errorbar(m3(indices),s3(indices),'-.k','linewidth',line_width);
    errorbar(m5(indices),s5(indices),'-c','linewidth',line_width);
    errorbar(m6(indices),s6(indices),'--m','linewidth',line_width);
    errorbar(m4(indices),s4(indices),'-g','linewidth',line_width);
    errorbar(m1(indices),s1(indices),'-r','linewidth',line_width);

    axis([1 15 0 1]);
    set(gca,'XTickLabel',{'3000','6000','9000','12000','15000'},'XTick',[3 6 9 12 15]);

    if my_cnt == 5
        xlabel('\bfTimestep','FontSize',12);
        legend('Degree','Betweenness','Closeness','Chen','EPN1','EPN2','Location','SouthEast');
    end
    if my_cnt == 1
        ylabel('\bfAverage of Opinions','FontSize',12);
    end
    % set(gca,'XScale','log');
    % title(['\bf(' char(my_ch) ')']);
else

    %% Real Networks
    disp('Real Networks');
    set(gcf, 'PaperPosition',[0.25 2.5 12 6]);
    figure(1);
    subplot(2,2,my_cnt);
    % runs = 1;
    line_width = 2;
    len = length(meanMajOpininos2Avg);
    % plot(1:length(meanMajOpininos2Avg),meanMajOpininos2Avg,'-k',1:length(meanMajOpininos3Avg),meanMajOpininos3Avg,'-.k',1:length(meanMajOpininos5Avg),meanMajOpininos5Avg,'-c',1:length(meanMajOpininos4Avg),meanMajOpininos4Avg,'--c',1:length(meanMajOpininos1Avg),meanMajOpininos1Avg,'-r','LineWidth',3);

    hold on;
    indices = 1:floor(len/COUNT):len;
    m2 = mean(meanMajOpininos2Avg);
    m3 = mean(meanMajOpininos3Avg);
    m5 = mean(meanMajOpininos5Avg);
    m4 = mean(meanMajOpininos4Avg);
    m1 = mean(meanMajOpininos1Avg);
    m6 = mean(meanMajOpininos6Avg);

    scale = 3;
    s2 = scale * std(meanMajOpininos2Avg);    %std(meanMajOpininos2Avg)/sqrt(runs);
    s3 = scale * std(meanMajOpininos3Avg);    %std(meanMajOpininos3Avg)/sqrt(runs);
    s5 = scale * std(meanMajOpininos5Avg);    %std(meanMajOpininos5Avg)/sqrt(runs);
    s4 = scale * std(meanMajOpininos4Avg);    %std(meanMajOpininos4Avg)/sqrt(runs);
    s1 = scale * std(meanMajOpininos1Avg);     %std(meanMajOpininos1Avg)/sqrt(runs);
    s6 = scale * std(meanMajOpininos6Avg);    %std(meanMajOpininos6Avg)/sqrt(runs);

    errorbar(m2(indices),s2(indices),'-k','linewidth',line_width);
    errorbar(m3(indices),s3(indices),'-.k','linewidth',line_width);
    errorbar(m5(indices),s5(indices),'-c','linewidth',line_width);
    errorbar(m6(indices),s6(indices),'--m','linewidth',line_width);
    errorbar(m4(indices),s4(indices),'-g','linewidth',line_width);
    errorbar(m1(indices),s1(indices),'-r','linewidth',line_width);

    axis([1 15 0 1]);
    set(gca,'XTickLabel',{'3000','6000','9000','12000','15000'},'XTick',[3 6 9 12 15]);

    if my_cnt == 2
        xlabel('\bfTimestep','FontSize',12);
        legend('Degree','Betweenness','Closeness','Chen','EPN1','EPN2','Location','SouthEast');
    end
    if my_cnt == 1
        ylabel('\bfAverage of Opinions','FontSize',12);
    end
    % set(gca,'XScale','log');
    % title(['\bf(' char(my_ch) ')']);

end
