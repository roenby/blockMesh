function b = pentagonInCircle(nth,nr2,r1,r2,z)

%b = pentagonInCircle(nth,nr2,r1,r2,z,prec) defines a cylindrical mesh
%where the outer rim is circular and the core is made up of a pentagon
%consisting of 5 quadrilaterals connected at the centre of the cylinder.
%There will be 10*nth cells around the perimeter of the cylinder.
%There will be nr2 radial points in the outer, non-pentagonal part of the
%cylinder. The cylinder will have radius r2 and the side length of the
%inner polygon will be r1. The vertical "floors" of the cylinder are
%defined in the array z.
%
%Johan Roenby, DHI Water & Environment

if nargin < 1
    close all
    nth = 3;
    nr2 = 6;
    r1 = 1;
    r2 = 2;
    z = [0 1];
end
prec = 1e-6;

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