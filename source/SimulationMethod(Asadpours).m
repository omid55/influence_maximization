%Omid55
%Simulation Method Function
function [ meanAllOpinions,meanAgOpininos,meanMajOpininos ] = SimulationMethod( x,N,InformedAgentsSize,MaximumSimulationSteps,mu,muPrime,XD,UI,alpha,sp,net,mode )

%% Approaches
if mode == 1
    %%% Random Selection Agents
    targets = [];
    nodesList = 1:N;
    nodesList(targets) = [];
    for i=1:InformedAgentsSize
        index = randi(length(nodesList));
        ia = nodesList(index);    % informed agents
        nodesList(index) = [];
        targets = [targets ia];
    end
else
    if mode == 2
        %%% Find High Degrees Agents
        degs = zeros(N,1);
        for i=1:N
            degs(i) = length(net(num2str(i)));
        end
        sortedDegs = sort(degs,'descend');
        H = InformedAgentsSize;    %H = 1;      % << CHECK HERE >>
        threshold = sortedDegs(H);
        indices = find(degs >= threshold);
        highDegs = indices(1:H)';
        targets = repmat(highDegs,1,InformedAgentsSize/H);
    else
        %%% Proposed Approach
        degs = zeros(N,1);
        for i=1:N
            degs(i) = length(net(num2str(i)));
        end

        % we want to find the nodes with small degree which have neighbors with big degrees
        deg_measure = zeros(N,1);
        for m=1:N
            if degs(m) > 1
                adjz = net(num2str(m));
                if degs(m) >= 2
                    deg_measure(m) = mean(degs(adjz)) / degs(m);
                end
            end
        end

        measure = (deg_measure./max(deg_measure));
        sortedMeasures = sort(measure,'descend');
        
        % Approach 1
        H = InformedAgentsSize;
        threshold = sortedMeasures(H);
        indices = find(measure >= threshold);
        highMeasures = indices(1:H)';
        targets = repmat(highMeasures,1,InformedAgentsSize/H);   
        
%         % Approach 2
%         targets = [];
%         list = 1:N;
%         while length(targets) < InformedAgentsSize
%             indices = find(measure == sortedMeasures(1));
%             sortedMeasures(1) = [];
%             goodIndx = indices(1);
%             measure(goodIndx) = [];
%             node = list(goodIndx);
%             list(goodIndx) = [];
%             targets = [targets ones(1,degs(node)+1)*node];    %targets = [targets ones(1,1)*node];
%         end
%         targets = targets(1:InformedAgentsSize);
        
    end
end

agents = N+1:N+InformedAgentsSize;
%x = 2 * rand(N+InformedAgentsSize,1) - 1;
x(agents) = XD;


%% Network Edge Adding
cnt = N+1;
for i=1:InformedAgentsSize
    nStr = num2str(targets(i));
    if net.isKey(nStr) == 0
        net(nStr) = cnt;
    else
        net(nStr) = [net(nStr) cnt];
    end
    net(num2str(cnt)) = targets(i);
    sp(targets(i),cnt) = 1;
    sp(cnt,targets(i)) = 1;
    cnt = cnt + 1;
end

% % --== if we can add informed agents to multiple nodes ==--
% cnt = 1;
% for i=1:length(targets)
%     nStr = num2str(targets(i));
%     leng = length(net(nStr));
%     for j=1:leng+1
%         if net.isKey(nStr) == 0
%             net(nStr) = N + cnt;
%         else
%             if any(net(nStr) == N + cnt)
%                 break;
%             end
%             net(nStr) = [net(nStr) N + cnt];
%         end
%         net(num2str(N + cnt)) = targets(i);
%         cnt = 1 + mod(cnt,InformedAgentsSize);
%     end
% end

%degree calculation
degs = zeros(N+InformedAgentsSize,1);
for i=1:N+InformedAgentsSize
    degs(i) = length(net(num2str(i)));
end


%% Simulation Body
allPos = [];
agPos = [];
majPos = [];
agentsOpinions = x(agents);
majorityOpinions = x(1:N);
meanAllOpinions = mean(x);
meanAgOpininos = mean(agentsOpinions);
meanMajOpininos = mean(majorityOpinions);

adjs = cell(N+InformedAgentsSize,1);
for i=1:N+InformedAgentsSize
    adjs{i} = net(num2str(i));
end

for l=1:MaximumSimulationSteps

%    %% Viewing The Graph
%     %%%
%     if mode == 3
%         labs=cell(1,InformedAgentsSize);
%         for ag=1:N+InformedAgentsSize
%             labs{ag} = [num2str(ag) ': ' num2str(x(ag))];
%         end
%         bg = biograph(sp,labs);
%         view(bg);
%         pause;
% 
%         % close all biograph windows
%         child_handles = allchild(0);
%         names = get(child_handles,'Name');
%         k = find(strncmp('Biograph Viewer', names, 15));
%         close(child_handles(k))
%     end
%     %%%
    
% % % % %     for i=1:N+InformedAgentsSize
% % % % %         
% % % % %         lastX = x;  
% % % % %         
% % % % %         adj = adjs{i};
% % % % %         if length(adj) == 0
% % % % %             continue;
% % % % %         end
% % % % %         
% % % % %         index = adj(randi(length(adj)));
% % % % %         Si = degs(i) ^ alpha;
% % % % %         Sj = degs(index) ^ alpha;
% % % % % 
% % % % %         Randomness = 0.1;       % << CHECK HERE FOR A SMALL VALUE FOR RANDOMNESS >>        Randomness = 0.0;       % << CHECK HERE FOR A SMALL VALUE FOR RANDOMNESS >>
% % % % %         if i > N    % if this is a informed agent            
% % % % %             if abs(lastX(i) - lastX(index)) <= UI 
% % % % %                 x(i) = StayInBound(lastX(i) + mu * (Sj/(Si+Sj)) * (lastX(index) - lastX(i)) + muPrime * (XD - lastX(i)));
% % % % %                 %x(i) = StayInBound(lastX(i) + mu * (lastX(index) - lastX(i)) + muPrime * (XD - lastX(i)));
% % % % %             end
% % % % %         else
% % % % %             uncertaintyI = rand(1) * Randomness + 1 - abs(lastX(i));
% % % % %             if abs(lastX(i) - lastX(index)) <= uncertaintyI
% % % % %                 x(i) = StayInBound(lastX(i) + mu * (Sj/(Si+Sj)) * (lastX(index) - lastX(i)));
% % % % %                 %x(i) = StayInBound(lastX(i) + mu * (lastX(index) - lastX(i)));
% % % % %             end
% % % % %         end
% % % % %         
% % % % %         if index > N            
% % % % %             if abs(lastX(i) - lastX(index)) <= UI
% % % % %                 x(index) = StayInBound(lastX(index) + mu * (Si/(Sj+Si)) * (lastX(i) - lastX(index)) + muPrime * (XD - lastX(index)));
% % % % %                 %x(index) = StayInBound(lastX(index) + mu * (lastX(i) - lastX(index)) + muPrime * (XD - lastX(index)));
% % % % %             end
% % % % %         else
% % % % %             uncertaintyIndex = rand(1) * Randomness + 1 - abs(lastX(index));
% % % % %             if abs(lastX(i) - lastX(index)) <= uncertaintyIndex
% % % % %                 x(index) = StayInBound(lastX(index) + mu * (Si/(Sj+Si)) * (lastX(i) - lastX(index)));
% % % % %                 %x(index) = StayInBound(lastX(index) + mu * (lastX(i) - lastX(index)));
% % % % %             end
% % % % %         end
% % % % %     end

   %% Asadpour's Way

    WMN = 0.3;
    WMS = 0.7;
    UM = 0.5;
    WIS = 0;
    WIN = 0.7;
    WID = 0.3;
    UI = 2;
    alpha = 0.25;
    XD = 1;    % Final Opinion
   
    for i=1:N+InformedAgentsSize
        lastX = x;

        adj = Adjacents(sp,i);
        if isempty(adj)
            continue;
        end

        index = 1 + floor(rand(1)*length(adj));

        if ~isempty(find(agents == i, 1))    %if this node is an informed agent
%             if abs(lastX(i) - lastX(index)) < UI
%                 r3 = 2 * rand(1) - 1;
%                 r4 = 2 * rand(1) - 1;
%                 r5 = 2 * rand(1) - 1;
%                 wiin = WIN * (1 + alpha * r3);
%                 wiis = WIS * (1 + alpha * r4);
%                 wiid = WID * (1 + alpha * r5);
%                 x(i) = wiin * lastX(index) + wiis * lastX(i) + wiid * XD;  
%             end
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

   %%

    positiveLeng = length(find(x>0));
    allPos = [allPos positiveLeng]; 
    agentsOpinions = x(agents);
    majorityOpinions = x(1:N);
    
    agPos = [agPos length(find(agentsOpinions>0))];
    majPos = [majPos length(find(majorityOpinions>0))];
    
    meanAllOpinions = [meanAllOpinions mean(x)];
    meanAgOpininos = [meanAgOpininos mean(agentsOpinions)];
    meanMajOpininos = [meanMajOpininos mean(majorityOpinions)];
    
%     if positiveLeng > 0.9*(N+InformedAgentsSize)
%         break;
%     end
end

end

