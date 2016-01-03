function A = ModularBA(m,n,k,k0,P_inter)

% ModularSW generates the adjacenecy matrix of modular scale-free networks 
% P_inter being the probability of inter-modular connections. 

A = zeros(m*n,m*n);
for i = 1:m;
    AA = BAnet(n, k, k0); 
    A((i-1)*n+1:i*n,(i-1)*n+1:i*n) = AA;
end;

for i = 1:m;
    for p = (i-1)*n+1:i*n-1;
        for q = p+1:i*n;
            R = rand; 
            if R < P_inter && A(p,q) == 1
                A(p,q) = 0;
                A(q,p) = 0;
                B = [p,q]; 
                B = B(randi(length(B)));
                BB = randi(m*n); 
                if B ~= BB
                    A(B,BB) = 1; 
                    A(B,BB) = 1; 
                end;
            end;
        end;
    end;
end;

end