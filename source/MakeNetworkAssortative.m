% Omid55
% we do not need to preserve scale free structure
function [ net ] = MakeNetworkAssortative( net , desiredAssortativity )

% %% You may delete this part
% N = 1000;
% averageDegree = 12;
% net = BAnet(N,averageDegree/2,averageDegree/2);
% %% You may delete this part

N = size(net,1);
it = 0;
cnt = 0;
finished = 0;

if desiredAssortativity > 0     % Assortative
    N = size(net,1);
    degs = full(sum(net));
    [~,idx] = sort(degs,'descend');
    half = floor(N/2);
    high = idx(1:half);
    low = idx(half+1:end);
    tmpNet = net;
    for i = 1 : length(high)
        for j = 1 : length(low)
            if tmpNet(high(i),low(j))
                tmpNet(high(i),low(j)) = 0;
                for k = 1 : length(high)
                    if tmpNet(high(i),high(k))==0 && i~=k
                        tmpNet(high(i),high(k)) = 1;
                        cnt = cnt + 1;
                    end
                end
            end
        end
    end
    cnt = cnt * desiredAssortativity;
    
    for i = 1 : length(high)
        for j = 1 : length(low)
            if net(high(i),low(j))
                net(high(i),low(j)) = 0;
                for k = 1 : length(high)
                    if net(high(i),high(k))==0 && i~=k
                        net(high(i),high(k)) = 1;
                        it = it + 1;
                        if it > cnt
                            ResultAssortativity = assortativity(net,0)
                            finished = 1;
                            break;
                        end
                    end
                end
                if finished == 1
                    break;
                end
            end
        end
        if finished == 1
            break;
        end
    end
else
    if desiredAssortativity < 0    % Disassortative
        degs = full(sum(net));
        [~,idx] = sort(degs,'descend');
        half = floor(N/2);
        high = idx(1:half);
        low = idx(half+1:end);
        tmpNet = net;
        for i = 1 : length(high)
            for j = 1 : length(high)
                if tmpNet(high(i),high(j))==1  && i~=j
                    tmpNet(high(i),high(j)) = 0;
                    for k = 1 : length(low)
                        if tmpNet(high(i),low(k)) == 0
                            tmpNet(high(i),low(k)) = 1;
                            cnt = cnt + 1;
                        end
                    end
                end
            end
        end
        cnt = cnt * -desiredAssortativity / 3.6;
        
        for i = 1 : length(high)
            for j = 1 : length(high)
                if net(high(i),high(j))==1  && i~=j
                    net(high(i),high(j)) = 0;
                    for k = 1 : length(low)
                        if net(high(i),low(k)) == 0
                            net(high(i),low(k)) = 1;
                            it = it + 1;
                            if it > cnt
                                ResultAssortativity = assortativity(net,0)
                                finished = 1;
                                break;
                            end
                        end
                    end
                    if finished == 1
                        break;
                    end
                end
            end
            if finished == 1
                break;
            end
        end
    else     % invalid type
        net = [];
    end
end

end

