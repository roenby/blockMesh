function [x,y] = rayPatch(L1,L2,s,a)

if nargin < 4
    a = 1;
end
    
if length(s) == 1 %Then s = number of cells
    s = cumsum([0, a.^[0:s-2]]);
end    
s = s/max(s);

for n = 1:size(L1,1)
    x(n,:) = lineSegment(L1(n,1),L2(n,1),s);
    y(n,:) = lineSegment(L1(n,2),L2(n,2),s);
end
x = x'; y = y';