%Omid55
function [  ] = RunAll(  )

%% Clear Everything
clc;
close all;
% addpath('/agbs/cluster/oaskaris/im/mb/matlab_bgl');


%% Run All
resPath = 'Results';
mkdir(resPath);

for netType = 1 : 1
    RunMain(resPath,netType);
end

disp('DONE');

end
