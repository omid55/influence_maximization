%Omid55
%ploting small sparse graphs with matlab
function [  ] = GraphDrawing( sp )

N = size(sp,1);
labs = cell(1,N);
for ag=1:N
    labs{ag} = num2str(ag);      % [num2str(ag) ': ' num2str(x(ag))];
end
bg = biograph(sp,labs);
view(bg);
pause;

% close all biograph windows
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
close(child_handles(k))

end

