%Omid55
%Forest Fire Network Creator Code
%gives N for number nodes and PF for Forward burning probability and PF for Backward burning probability.
function [ sp ] = ForestFireCreator( N,PF,PB )

initCliqueSize = 3;
sp = sparse(ones(initCliqueSize) - diag(ones(initCliqueSize,1)));   % clique creator

for i=initCliqueSize+1:N
    ambassador = 1 + floor( rand(1) * (i-1) );
    sp(i,ambassador) = 1;
    sp(ambassador,i) = 1;
    nodes = getRecursivelyNodes(ambassador,sp,PF,PB);
    nodes(find(nodes == i))=[];
    for j=1:length(nodes)
        sp(i,nodes(j)) = 1;
        sp(nodes(j),i) = 1;
    end
end

end

