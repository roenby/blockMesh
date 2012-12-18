function b = smoothedPentagonBlock(p1,p2,p3,p4,p5,n,z,prec)

%b = smoothedPentagonBlock(p1,p2,p3,p4,p5,n,z,prec) generates a pentagonal
%block with corners at the points p1-p5. There will be 2*n cells along each
%side of the pentagon. The vertical "floors" of the pentagon block are
%specified in the array z. Two points are regarded as equal if they are
%within a distance of prec to each other.
%
%Johan Roenby, DHI Water & Environment

%Center point in pentagon
p0 = faceCentre([p1; p2; p3; p4; p5]);
%Middle points of sides of pentagon
p12 = 1/2*(p1 + p2);
p23 = 1/2*(p2 + p3);
p34 = 1/2*(p3 + p4);
p45 = 1/2*(p4 + p5);
p15 = 1/2*(p1 + p5);

%Blocks
[x,y] = smoothedPatch(p1,p12,p0,p15,n,n,prec);
b = repeat2DMesh(x,y,z);
[x,y] = smoothedPatch(p2,p23,p0,p12,n,n,prec);
b2 = repeat2DMesh(x,y,z);
b = mergeBlocks(b,b2,prec);
[x,y] = smoothedPatch(p3,p34,p0,p23,n,n,prec);
b2 = repeat2DMesh(x,y,z);
b = mergeBlocks(b,b2,prec);
[x,y] = smoothedPatch(p4,p45,p0,p34,n,n,prec);
b2 = repeat2DMesh(x,y,z);
b = mergeBlocks(b,b2,prec);
[x,y] = smoothedPatch(p5,p15,p0,p45,n,n,prec);
b2 = repeat2DMesh(x,y,z);
b = mergeBlocks(b,b2,prec);

function p = faceCentre(p)

if size(p,2) == 2;
    p = [p, zeros(size(p,1),1)];
end

%Initial guess of face centres
x = sum(p,1)/size(p,1);

%Calculating surface normal and face centre by triangulation with
%initially guessed face centre
%Sf = zeros(size(f,1),3);
Cf = zeros(1,3);
sumAbsSf = 0;
for k = 1:size(p,1)
    x1 = p(k,:);
    kPlus1 = mod(k,size(p,1))+1;
    x2 = p(kPlus1,:);
    Sk = 1/2*cross(x2-x1,x-x1);
%    Sf = Sf + Sk;
    Cf = Cf + 1/3*norm(Sk).*(x1 + x2 + x);
    sumAbsSf = sumAbsSf + norm(Sk); 
end
p = Cf(1:2)./sumAbsSf;