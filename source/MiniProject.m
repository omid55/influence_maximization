%Omid55
%Complex Dynamical Networks' Mini Project
function [  ] = MiniProject(  )


%% Initiallization
populationSize = 1000;    
informedAgentsSize = 100; 
averageDegree = 30;
MaximumSimulationSteps = 100000;

WMN = 0.3;
WMS = 0.7;
UM = 0.5;
WIS = 0;
WIN = 0.7;
WID = 0.3;
UI = 2;
alpha = 0.25;
XD = 1;    % Final Opinion


%% Creating Network (Erdos Renyi)
net = BarabasiGraphCreator(N,4,4);


%% Calculating Values
agents = [];
nodesList = 1:populationSize;
for i=1:informedAgentsSize
    index = 1 + floor(rand(1)*length(nodesList));
    ia = nodesList(index);    % informed agents
    nodesList(index) = [];
    agents = [agents ia];
end

x = [];
for i=1:populationSize
    x(i) = 2*rand(1)-1;
end


% you can comment it for another senario
for i=1:informedAgentsSize
    x(agents(i)) = XD;
end
% you can comment it for another senario



%% Simulation Body
allPos = [];
agPos = [];
majPos = [];
meanAllOpinions = [];
meanAgOpininos = [];
meanMajOpininos = [];

for l=1:MaximumSimulationSteps

    for i=1:populationSize
        lastX = x;

        adj = Adjacents(sp,i);
        if length(adj) == 0
            continue;
        end

        index = 1 + floor(rand(1)*length(adj));

        if length(find(agents == i)) ~= 0    %if this node is an informed agent
            if abs(lastX(i) - lastX(index)) < UI
                r3 = 2 * rand(1) - 1;
                r4 = 2 * rand(1) - 1;
                r5 = 2 * rand(1) - 1;
                wiin = WIN * (1 + alpha * r3);
                wiis = WIS * (1 + alpha * r4);
                wiid = WID * (1 + alpha * r5);
                x(i) = wiin * lastX(index) + wiis * lastX(i) + wiid * XD;  
            end
        else    % if this is a majority
            if abs(lastX(i) - lastX(index)) < UM
                r1 = 2 * rand(1) - 1;
                r2 = 2 * rand(1) - 1;
                wimn = WMN * (1 + alpha * r1);
                wims = WMS * (1 + alpha * r2);
                x(i) = wimn * lastX(index) + wims * lastX(i);
            end
        end
    end

    positiveLeng = length(find(x>0));
    
    allPos = [allPos positiveLeng]; 
    
    agentsOpinions = x(agents);
    maj = 1:populationSize;
    maj = setdiff(maj,agents);
    majorityOpinions = x(maj);
    
    agPos = [agPos length(find(agentsOpinions>0))];
    majPos = [majPos length(find(majorityOpinions>0))];
    
    meanAllOpinions = [meanAllOpinions mean(x)];
    meanAgOpininos = [meanAgOpininos mean(agentsOpinions)];
    meanMajOpininos = [meanMajOpininos mean(majorityOpinions)];
    
    if positiveLeng > 0.9*populationSize
        break;
    end
end

figure,plo = plot(allPos);
xlabel('Timestep');
ylabel('Population with Positive Opinion');
title('All Society Temporal Evolution of the Opinion Formation Process');
legend(plo,'Whole Society');

figure,p1 = plot(agPos,'r');
xlabel('Timestep');
ylabel('Informed Agents & Majority Population with Positive Opinion');
title('Temporal Evolution of the Opinion Formation Process');
hold on;
p2 = plot(majPos,'b');
hold on;
legend([p1,p2],'Informed Agents','Majority');



figure,plo = plot(meanAllOpinions);
xlabel('Timestep');
ylabel('Average of Opinions');
title('All Society Temporal Evolution of the Opinion Formation Process');
legend(plo,'Whole Society');

figure,p1 = plot(meanAgOpininos,'r');
xlabel('Timestep');
ylabel('Average of Opinions');
title('Temporal Evolution of the Opinion Formation Process');
hold on;
p2 = plot(meanMajOpininos,'b');
hold on;
legend([p1,p2],'Informed Agents','Majority');

end

