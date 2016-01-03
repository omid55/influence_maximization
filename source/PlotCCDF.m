figure;
plot(sort(d0), linspace(1-1/length(d0), 0, length(d0)),'-r','LineWidth',2);
hold on;
plot(sort(d1), linspace(1-1/length(d1), 0, length(d1)),'-b','LineWidth',2);
hold on;
plot(sort(d2), linspace(1-1/length(d2), 0, length(d2)),'-g','LineWidth',2);
hold on;
plot(sort(d3), linspace(1-1/length(d3), 0, length(d3)),'-c','LineWidth',2);
legend('Degree','Betweenness','Proposed1','Proposed2');
ylabel('\bfP(X>o)');
xlabel('\bfOpinion');

figure;
bin = -0.2:0.01:0.2;
x = (bin(1:length(bin)-1)+bin(2:length(bin)))/2;
plot(x,MyBinsCounting(d0,bin),x,MyBinsCounting(d1,bin),x,MyBinsCounting(d2,bin),x,MyBinsCounting(d3,bin),'LineWidth',2);
legend('Degree','Betweenness','Proposed1','Proposed2');
xlabel('\bfOpinion');
ylabel('\bfVolume');


% figure;
% plot(sort(d0(end,:)), linspace(1-1/length(d0(end,:)), 0, length(d0(end,:))),'-r');
% hold on;
% plot(sort(d1(end,:)), linspace(1-1/length(d1(end,:)), 0, length(d1(end,:))),'-b');
% hold on;
% plot(sort(d2(end,:)), linspace(1-1/length(d2(end,:)), 0, length(d2(end,:))),'-g');
% hold on;
% plot(sort(d3(end,:)), linspace(1-1/length(d3(end,:)), 0, length(d3(end,:))),'-c');
% legend('Degree','Betweenness','Proposed1','Proposed2');
% ylabel('\bfP(X>o)');
% xlabel('\bfOpinion');
% 
% 
% figure;
% plot(1:size(d0,1),mean(d0'),'-o',1:size(d1,1),mean(d1'),'-p',1:size(d2,1),mean(d2'),'-x',1:size(d3,1),mean(d3'),'-^');
% legend('Degree','Betweenness','Proposed1','Proposed2');
% 
% 
% figure;
% network_file = 'BA100.txt';
% a = importdata(network_file);
% MM = max(max(a(:,1)),max(a(:,2)));
% Graph = sparse(a(:,1),a(:,2),1,MM,MM);
% plot(mean(d3(:,find(Graph(77,:)==1))'));