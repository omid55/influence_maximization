%Omid55
%Simulation Method Function
function [ meanMajOpinions,lastOpinions,fullDistribution,followers ] = SimulationMethod( x,InformedAgentsSize,MaximumSimulationSteps,mu,A,mode,Randomness,alpha,beta,betweennesses,closenesses )

%% Modes
N = size(x,1);
XD = 1;
indegs = sum(A);
outdegs = sum(A');
%alpha = 0;
majorityMeanOpinions = mean(x);
initFol = length(find(x>0)) / length(x);

if mode == 1
%     %% Random Selection Agents
%     targets = datasample(1:N,InformedAgentsSize,'Replace',false);
    
   %% EPN2
    % we want to find the nodes with small degree which have neighbors with large degrees
     EPN2_measure = zeros(N,1);
     eps = 0.0001;
     for i = 1 : N
         f1 = ((outdegs(i) + eps) / (indegs(i) + eps));
         ss = 0;
         for j = 1 : outdegs(i)
            ss =  ss + ( beta * indegs(j) + (1-beta) * mean(indegs(A(i,:) == 1)) );
         end
         ss = ss / outdegs(i);
         f2 = ss / max(indegs);
         EPN2_measure(i) = f1 * f2;
     end
     [~,EP2] = sort(EPN2_measure,'descend');
     targets = EP2(1:InformedAgentsSize)';
else
    if mode == 2
        %%% Find High Degrees Agents
        [~,idxx] = sort(indegs,'descend');
        targets = idxx(1:InformedAgentsSize)';
    else
        if mode == 3
            [~,idddx] = sort(betweennesses,'descend');
            targets = idddx(1:InformedAgentsSize)';
%             % Inverse of degree
%             [~,idxx] = sort(indegs,'ascend');
%             targets = idxx(1:InformedAgentsSize)';
        else
            if mode == 4
           %% EPN1
                EPN1_measure = zeros(N,1);
                eps = 0.0001;
                for i = 1 : N
                    f1 = ((outdegs(i) + eps) / (indegs(i) + eps));
                    f2 = (mean(indegs(A(i,:) == 1)) / max(indegs));
                    EPN1_measure(i) = f1 * f2;
                end
                [~,EP1] = sort(EPN1_measure,'descend');
                targets = EP1(1:InformedAgentsSize)';
            else
                if mode == 5
                    [~,xx] = sort(closenesses,'descend');
                    targets = xx(1:InformedAgentsSize)';
                else
                    % Chen measure
                    if mode == 6
                        NN = zeros(N,1);
                        for i = 1 : N
                            list = find(A(:,i) == 1);
                            list2 = [];
                            for j = 1 : length(list)
                                list2 = union(list2,find(A(:,list(j))==1));
                            end
                            NN(i) = length(setdiff(union(list,list2),i));
                        end
                        chen = zeros(N,1);
                        for i = 1 : N
                            addj = find(A(:,i) == 1);
                            for j = 1 : length(addj)
                                chen(i) = chen(i) + sum(NN(A(:,addj(j)) == 1));
                            end
                        end
                        
                        [~,xxx] = sort(chen,'descend');
                        targets = xxx(1:InformedAgentsSize)';
                    end
                end
            end
        end
    end
end

agents = N+1:N+InformedAgentsSize;
x(agents) = XD;
% AGENTS_ADDED = 1


%% Network Edge Adding
for i=1:length(targets)
    A(targets(i),N+i) = 1;
    A(N+i,targets(i)) = 1;
end
net = containers.Map;
for i = 1 : N+InformedAgentsSize
    net(num2str(i)) = find(A(i,:) == 1);
end
% indegs = sum(A);


%% Simulation Body
[meanMajOpinions,lastOpinions,fullDistribution,followers] = MySimulationBodyMethod(x,InformedAgentsSize,MaximumSimulationSteps,mu,A,Randomness,alpha);
meanMajOpinions = [majorityMeanOpinions meanMajOpinions];
followers = [initFol followers];


%% Simulation Body
% % Asadpour's Way
% WMN = 0.3;
% WMS = 0.7;
% %UM = 0.5;
% WIS = 0;
% WIN = 0.7;
% WID = 0.3;
% UI = 2;
% alpha = 0.25;
% XD = 1;    % Final Opinion
% meanMajOpininos = zeros(MaximumSimulationSteps+1,size(x,2));
% meanMajOpininos(1,:) = mean(x(1:N));
% for l=1:MaximumSimulationSteps
%     
%     lastX = x;
%     for i=1:N+InformedAgentsSize
%         
%         adj = net(num2str(i));
%         if isempty(adj)uj
%             continue;
%         end
% 
%         idx = randi(length(adj));
%         index = adj(idx);
% 
%         if i > N
%             if abs(lastX(i) - lastX(index)) < UI
%                 r3 = 2 * rand(1) - 1;
%                 r4 = 2 * rand(1) - 1;
%                 r5 = 2 * rand(1) - 1;
%                 wiin = WIN * (1 + alpha * r3);
%                 wiis = WIS * (1 + alpha * r4);
%                 wiid = WID * (1 + alpha * r5);
%                 x(i) = wiin * lastX(index) + wiis * lastX(i) + wiid * XD;  
%             end
%         else    % if this is a majority
%             uncertaintyI = rand(1) * Randomness + 1 - abs(lastX(i));
%             if abs(lastX(i) - lastX(index)) < uncertaintyI
%                 r1 = 2 * rand(1) - 1;
%                 r2 = 2 * rand(1) - 1;
%                 wimn = WMN * (1 + alpha * r1);
%                 wims = WMS * (1 + alpha * r2);
%                 x(i) = wimn * lastX(index) + wims * lastX(i);
%             end
%         end
%     end
% 
%     meanMajOpininos(l+1,:) = mean(x(1:N));
% end

% % My way
% indegs = sum(A);
% net = containers.Map;
% for i = 1 : N
%     net(num2str(i)) = find(A(i,:) == 1);
% end
% meanMajOpininos = zeros(MaximumSimulationSteps+1,size(x,2));
% meanMajOpininos(1,:) = mean(x(1:N));
% for l=1:MaximumSimulationSteps
%     
%     lastX = x;
%     
%     for i = 1 : N
%         adj = net(num2str(i));
%         uncertaintyI = rand(1) * Randomness + 1 - abs(lastX(i));
%         op = 0;
%         SPs = 0;
%         for k = 1 : length(adj)
%             j = adj(k);
%             SPj = indegs(j) ^ alpha;
%             if j > N
%                 SPj = 100;
%                 if rand() < uncertaintyI
%                     newOpinionForAgent = lastX(i) + uncertaintyI;
%                     if newOpinionForAgent > 1
%                         newOpinionForAgent = 1;
%                     end
%                     op = op + SPj * newOpinionForAgent;
%                     SPs = SPs + SPj;
%                 end
%             else
%                 if abs(lastX(i) - lastX(j)) <= uncertaintyI
%                     op = op + SPj * lastX(j);
%                     SPs = SPs + SPj;
%                 end
%             end
% %             if abs(lastX(i) - lastX(j)) <= uncertaintyI
% %                 op = op + SPj * lastX(j);
% %                 SPs = SPs + SPj;
% %             end
%         end
%         if SPs == 0
%             continue;
%         end
%         x(i) = op / SPs;
%     end
% 
%     meanMajOpininos(l+1,:) = mean(x(1:N));
%     
% 
% end


% %Based on RA
% meanMajOpinions = zeros(MaximumSimulationSteps+1,size(x,2));
% meanMajOpinions(1,:) = mean(x(1:N));
% us = zeros(N,1);
% for i = 1 : N
%     us(i) = 1 / (1+exp(indegs(i)));
% end
% for it=1:MaximumSimulationSteps
%     lastX = x;
%     for i = 1 : N
%         adj = net(num2str(i));
%         ui = (max(indegs) - indegs(i)) / max(indegs);
%         op = 0;
%         we = 0;
%         for k = 1 : length(adj)
%             j = adj(k);
%             uj = (max(indegs) - indegs(j)) / max(indegs);
%             if j > N
%                 uj = 10^-10;
%             end
%             hij = min(lastX(i)+ui,lastX(j)+uj) - max(lastX(i)-ui,lastX(j)-uj);
%             if hij > uj
%                 op = op + (hij/uj - 1) * (lastX(j) - lastX(i));
%                 we = we + (hij/uj - 1);
%             end
%         end
%         if we == 0
%             continue;
%         end
%         x(i) = op / we;
%     end
%     meanMajOpinions(it+1,:) = mean(x(1:N));
% end


% %Good approach
% meanMajOpinions = zeros(MaximumSimulationSteps+1,size(x,2));
% meanMajOpinions(1,:) = mean(x(1:N));
% for it=1:MaximumSimulationSteps
%     lastX = x;
%     for i = 1 : N
%         adj = net(num2str(i));
%         SPi = indegs(i) ^ alpha;
%         ui = (max(indegs) - indegs(i)) / max(indegs) + Randomness * rand(1);
%         op = 0;
%         SPs = 0;
%         for k = 1 : length(adj)
%             j = adj(k);
%             if j > N && rand() < ui
%                 SPj = N * 10;
%                 op = op + SPj * lastX(j);
%                 SPs = SPs + SPj;
%             else
%                 if abs(lastX(i) - lastX(j)) < ui
%                     SPj = indegs(j) ^ alpha;
%                     op = op + SPj * lastX(j);
%                     SPs = SPs + SPj;
%                 end
%             end
%         end
%         if SPs == 0
%             continue;
%         end
%         SPs = SPs + SPi;
%         op = op + SPi * lastX(i);
%         x(i) = op / SPs;
%     end
%     meanMajOpinions(it+1,:) = mean(x(1:N));
% end


end


