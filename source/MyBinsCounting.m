% Omid55
% x will be the median of each ranges and y will be the count
function [ y ] = MyBinsCounting ( times , binranges )

[bincounts] = histc(times,binranges);
y = bincounts(1:end-1);
y = y / max(y);

end

