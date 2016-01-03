function A = ModularWS(m,n,k,P_intra,P_inter)

% ModularWS generates the adjacenecy matrix of modular small-world networks 
% P_inter being the ptobability of inter-modular connections. 

A = zeros(m*n,m*n);
for i = 1:m;
    AA = WattsStrogatzCreator(n, k, P_intra); 
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