function L = lineSegment(p1,p2,s)

% If length(s) > 1 it is interpreted as the points from p1 to p2.
% Else s = number of bins, so number of points will be one more

p1 = p1(:)';
p2 = p2(:)';

if length(s) == 1
%    s = (0:s-1)'/(s-1);
    s = (0:s)'/s;
else
    s = s(:);
end

L = ones(size(s))*p1 + s*(p2-p1);