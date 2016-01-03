%Omid55
function [  ] = ResultGenerator ( ResultsPath,ttype )
clc;
close all;


%% Artificial Figures
if ttype == 1
    
    disp('Art');
    ResultsPath
    COUNT = 20;
    fig = figure;

    my_cnt = 1;
    my_ch = 'a';
    %ResultsPath = 'ServerResultsAlpha0'   %'ServerResultsAlpha1';

    load([ResultsPath '/BA/FullData.mat']);
    %load([ResultsPath '\BA\FullData.mat']);
    resgen;
    title(['\bf(' char(my_ch) '_' num2str(my_cnt) ')']);
    my_cnt = my_cnt + 1;

    load([ResultsPath '/ER/FullData.mat']);
    % load([ResultsPath '\ER\FullData.mat']);
    resgen;
    title(['\bf(' char(my_ch) '_' num2str(my_cnt) ')']);
    my_cnt = my_cnt + 1;

    load([ResultsPath '/ModularBA/FullData.mat']);
    % load([ResultsPath '\ModularBA\FullData.mat']);
    resgen;
    title(['\bf(' char(my_ch) '_' num2str(my_cnt) ')']);
    my_cnt = my_cnt + 1;

    load([ResultsPath '/WS/FullData.mat']);
    % load([ResultsPath '\WS\FullData.mat']);
    resgen;
    title(['\bf(' char(my_ch) '_' num2str(my_cnt) ')']);
    my_cnt = my_cnt + 1;

    load([ResultsPath '/ModularWS/FullData.mat']);
    % load([ResultsPath '\ModularWS\FullData.mat']);
    resgen;
    title(['\bf(' char(my_ch) '_' num2str(my_cnt) ')']);
    my_cnt = my_cnt + 1;

    load([ResultsPath '/FF/FullData.mat']);
    % load([ResultsPath '\FF\FullData.mat']);
    resgen;
    title(['\bf(' char(my_ch) '_' num2str(my_cnt) ')']);

    saveas(fig, [ResultsPath '.fig']);
    % print('-dpng', '-r1200', [ResultsPath '.png']);

else
    
    %% Real Figures
    disp('Real');
    ResultsPath
    COUNT = 20;
    fig = figure;

    my_cnt = 1;
    my_ch = 'a';

    load([ResultsPath '/Advogato/FullData.mat']);
    resgen;
    title(['\bf(' char(my_ch) '_' num2str(my_cnt) ')']);
    my_cnt = my_cnt + 1;

    load([ResultsPath '/Hamster/FullData.mat']);
    resgen;
    title(['\bf(' char(my_ch) '_' num2str(my_cnt) ')']);
    my_cnt = my_cnt + 1;

    load([ResultsPath '/HamsterFriendships/FullData.mat']);
    resgen;
    title(['\bf(' char(my_ch) '_' num2str(my_cnt) ')']);
    my_cnt = my_cnt + 1;

%     load([ResultsPath '/Facebook/FullData.mat']);
    load([ResultsPath '/Elec/FullData.mat']);
    resgen;
    title(['\bf(' char(my_ch) '_' num2str(my_cnt) ')']);
    my_cnt = my_cnt + 1;

    saveas(fig, [ResultsPath '.fig']);

end


end
