%Omid55
%Simulation Method Function
function [ meanMajOpininos ] = SimulationMethod( x,InformedAgentsSize,MaximumSimulationSteps,mu,A,mode,Randomness )

%% Modes
N = size(x,1);
%majorityMeanOpinions = mean(x);
XD = 1;
degs = sum(A);
alpha = 2;

if mode == 1
    %%% Random Selection Agents
    targets = datasample(1:N,InformedAgentsSize,'Replace',false);
else
    if mode == 2
        %%% Find High Degrees Agents
        [~,idxx] = sort(degs,'descend');
        targets = idxx(1:InformedAgentsSize)';
    else
        if mode == 3
            % Inverse of degree
            [~,idxx] = sort(degs,'ascend');
            targets = idxx(1:InformedAgentsSize)';
        else
            %%% Proposed Approach

            % we want to find the nodes with small degree which have neighbors with big degrees
            deg_measure = zeros(N,1);
            for m=1:N
                if degs(m) >= 2
                    adjz = A(m,:)==1;
                    deg_measure(m) = mean(degs(adjz)) / degs(m);
                end
            end
            measure = deg_measure./max(deg_measure);
            [~,idddx] = sort(measure,'descend');

%         % with betweenness
%         measure = betweenness_centrality(sparse(A));
%         sortedMeasures = sort(measure,'descend');
        
        %switch approach
            
            %case 1
            % Approach 1
            targets = idddx(1:InformedAgentsSize)';  
        
%             %case 2
%             % Approach 2
%             targets = [];
%             list = 1:N;
%             while length(targets) < InformedAgentsSize
%                 indices = find(measure == sortedMeasures(1));
%                 sortedMeasures(1) = [];
%                 goodIndx = indices(1);
%                 measure(goodIndx) = [];
%                 node = list(goodIndx);
%                 list(goodIndx) = [];
%                 targets = [targets ones(1,degs(node)+1)*node];
%             end
%             targets = targets(1:InformedAgentsSize);
            
        end
    end
end

agents = N+1:N+InformedAgentsSize;
%x = 2 * rand(N+InformedAgentsSize,1) - 1;
x(agents) = XD;


%% Network Edge Adding
for i=1:length(targets)
    A(targets(i),N+i) = 1;
    A(N+i,targets(i)) = 1;
end
% A(targets,N+1:N+InformedAgentsSize) = 1;    % this is NOT TRUE    << CHECK HERE >>
% A(N+1:N+InformedAgentsSize,targets) = 1;     % this is NOT TRUE     << CHECK HERE >>


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

% %degree calculation
% degs = zeros(N+InformedAgentsSize,1);
% for i=1:N+InformedAgentsSize
%     degs(i) = length(net(num2str(i)));
% end


%% Simulation Body
% meanMajOpininos = MySimulationBodyMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,Randomness);
% meanMajOpininos = [majorityMeanOpinions meanMajOpininos];

degs = sum(A);
meanMajOpininos = zeros(MaximumSimulationSteps+1,size(x,2));
meanMajOpininos(1,:) = mean(x(1:N));

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

    %%%i=randi(N);
    for i=1:N
        lastX = x;  
        
        adj = find(A(i,:)==1);
        if length(adj) == 0
            continue;
        end
                
%         %%
%         Randomness = 0.0;
%         removeList = [];
%         uncertaintyI = rand(1) * Randomness + 1 - abs(lastX(i));
%         for j=1:length(adj)
%             if abs(lastX(i) - lastX(adj(j))) > uncertaintyI
%                 removeList = [removeList j];
%             end
%         end
%         adj(removeList) = [];
%         if length(adj) == 0
%             continue;
%         end
% 
%         diff = 0;
%         for j=1:length(adj)
%             index = adj(j);
%             Si = degs(i) ^ alpha;
%             Sj = degs(index) ^ alpha;
%             diff = diff + mu * (Sj/(Sj+Si)) * (lastX(index) - lastX(i));
%         end
%         x(i) = StayInBound(lastX(i) + diff/length(adj));
%         %%
        
        index = adj(randi(length(adj)));
        Si = degs(i) ^ alpha;
        Sj = degs(index) ^ alpha;
        
        if index > N
            Sj = 50;
        end

        %Randomness = 0.01;       % << CHECK HERE FOR A SMALL VALUE FOR RANDOMNESS >>        Randomness = 0.0;       % << CHECK HERE FOR A SMALL VALUE FOR RANDOMNESS >>        
        uncertaintyI = rand(1) * Randomness + 1 - abs(lastX(i));
        if abs(lastX(i) - lastX(index)) <= uncertaintyI
           x(i) = StayInBound(lastX(i) + mu * (Sj/(Si+Sj)) * (lastX(index) - lastX(i)));
           %x(i) = StayInBound(lastX(i) + mu * (lastX(index) - lastX(i))); 
        end
        
        if index <= N
            uncertaintyIndex = rand(1) * Randomness + 1 - abs(lastX(index));
            if abs(lastX(i) - lastX(index)) <= uncertaintyIndex
                x(index) = StayInBound(lastX(index) + mu * (Si/(Sj+Si)) * (lastX(i) - lastX(index)));
                %x(index) = StayInBound(lastX(index) + mu * (lastX(i) - lastX(index)));
            end
        end
    end
    
%     meanAllOpinions(l+1,:) = mean(x);
%     meanAgOpininos(l+1,:) = mean(x(agents));
    meanMajOpininos(l+1,:) = mean(x(1:N));
    
%     if positiveLeng > 0.9*(N+InformedAgentsSize)
%         break;
%     end
end

end

