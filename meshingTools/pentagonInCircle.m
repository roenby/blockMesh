function b = pentagonInCircle(nth,nr2,r1,r2,z,prec)

if nargin < 1
    close all
    nth = 3;
    nr2 = 6;
    r1 = 1;
    r2 = 2;
    z = [0 1];
    prec = 1e-6;
end

pc = r1*exp(2*pi*1i*(0:4)/5);
p = zeros(length(pc),2);
for n = 1:length(pc)
    p(n,:) = [real(pc(n)) imag(pc(n))];
end
b = smoothedPentagonBlock(p(1,:),p(2,:),p(3,:),p(4,:),p(5,:),nth,z,prec);

L0 = lineSegment(p(1,:),p(2,:),2*nth);
th = 2*pi/5*(0:2*nth)/(2*nth);
L1 = r2*[cos(th(:)) sin(th(:))];
[x,y] = rayPatch(L0,L1,nr2);
b2 = repeat2DMesh(x,y,z);
for n = 0:4
    b2 = rotBlock(b2,[0 0 2*pi/5]);
    b = mergeBlocks(b,b2,prec);
end

if nargin < 1
    plot3(b.points(:,1),b.points(:,2),b.points(:,3),'.')
    view(2)
    axis equal
    plotMesh(b)
end