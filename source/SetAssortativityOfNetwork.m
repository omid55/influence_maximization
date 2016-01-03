% Omid55
% Based on (Milo et al., 2003; Payne and Eppstein, 2009)

% A standard edge-swapping method was employed for tuning the assortativity of a network to a desired value.
% In each iteration of this method two edges i->j and x->y were selected at random. Edges were then swapped,
% resulting in two new edges i->y and x->j. These new edges replaced i->j and x->y if they changed assortativity
% in the desired direction. Otherwise, the new edges were discarded and the old edges were kept. Such edge 
% swaps preserved the in- and out-degrees of all nodes involved, thereby keeping the degree distribution intact.

function [ net ] = SetAssortativityOfNetwork( net , desiredAssortativity )

%% Init
addpath('BCT');
isDirectedNetwork = 1;
eps = 0.1;
N = size(net,1);

%% Calculation
currentAssortativity = assortativity(net,isDirectedNetwork);
distance = abs(currentAssortativity - desiredAssortativity);
while distance > eps
    indices = randsample(N,2);
    i = indices(1);
    x = indices(2);
    Ni = find(net(i,:) == 1);
    Nx = find(net(x,:) == 1);
    if isempty(Ni) || isempty(Nx)
        continue;
    end
    j = Ni(randi(length(Ni),1));
    y = Nx(randi(length(Nx),1));
    net(i,j) = 0;
    net(x,y) = 0;
    net(i,y) = 1;
    net(x,j) = 1;
    newDistance = abs(assortativity(net,isDirectedNetwork) - desiredAssortativity);
    if newDistance < distance
        distance = newDistance;
    else
        net(i,j) = 1;
        net(x,y) = 1;
        net(i,y) = 0;
        net(x,j) = 0;
    end
end

end

