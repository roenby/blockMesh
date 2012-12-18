function b = rescaleBlock(b,s)

if length(s) == 1
    s = s*ones(size(b.points,2));
end

for n = 1:size(b.points,2)
    b.points(:,n) = s(n)*b.points(:,n);
end