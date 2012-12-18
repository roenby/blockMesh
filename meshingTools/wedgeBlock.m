function b = wedgeBlock(Nx,Ny,nz)

%b = wedgeBlock(Nx,Ny) generates an Nx by Ny by nz*Nx block which is 
%triangular in the xz-plane so that the mesh goes from [0:Nx] in the x
%direction follows the line z = nz*x in and goes from [0:nz*Nx] in the z
%direction. In the y direction it goes from 0 to Ny.
%Run the code without any arguments to see a plot.
%
%Johan Roenby, DHI Water & Environment

if nargin < 1
    Nx = 5;
    Ny = 2;
end
if nargin < 3
    nz = 1;
end
prec = 1e-6;

b = block(Nx,Ny,nz);
b.points(:,3) = b.points(:,3);
b.points(:,1) = b.points(:,1);
for n = 1:Nx-1
    b2 = block(Nx-n,Ny,nz);
    b2.points(:,3) = b2.points(:,3) + n*nz;
    b2.points(:,1) = b2.points(:,1) + n;
    b = mergeBlocks(b,b2,prec);
end

P = b.points;
x = [0 0 0];
n = [-nz 0 1];
P = moveToPlane(P,x,n);
b.points = P;

if nargin < 1
    clf; plotMesh(b)
end

function P = moveToPlane(P,x,n)

X = repmat(x(:)',size(P,1),1);
N = repmat(n(:)'/norm(n),size(P,1),1);
para = dot(P-X,N,2);
ind = find(para > 0);
perp = cross(N(ind,:),cross(P(ind,:)-X(ind,:),N(ind,:)));
%alpha = (P(ind,3) - X(ind,3))./perp(:,3);
%Pnew = repmat(alpha,1,3).*perp + X(ind,:);
%alpha(Pnew(:,1)-P(ind,1) >= 1) = 1;
%P(ind,:) = repmat(alpha,1,3).*perp + X(ind,:);
P(ind,:) = perp + X(ind,:);