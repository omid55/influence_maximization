% data = [1 2; 1 3;1 4; 1 5; 1 6; 1 7; 1 8; 1 9; 6 10; 10 11; 10 23; 11 12; 11 21; 11 23; 12 14; 12 13; 13 14; 13 15; 13 22; 14 15; 14 23; 15 16; 16 17; 16 18; 16 22; 17 18; 17 20; 17 21; 18 19; 18 22; 19 20; 19 21; 20 21; 20 23; 22 23; 2 3; 3 4; 7 8; 8 9];
% sp = sparse(data(:,1),data(:,2),1,23,23);
% A = full(sp);
% A = max(A,A');


% CHEN
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
