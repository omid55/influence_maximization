%Omid55
%ER
function [ net,sp ] = ErdosRenyi( N,ErdosRenyiProbablity )

mat = rand(N,N);
index = find(mat < ErdosRenyiProbablity);
A = zeros(N,N);
A(index) = 1;
sp = sparse(A);
sp = max(sp,sp');
net = CreateMap(sp);

end

