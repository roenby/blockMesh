function b = wedgeBlock(Nx,Ny,prec)

if nargin < 1
    Nx = 5;
    Ny = 3;
    prec = 1e-6;
end

nz = 0;
b = block(Nx,Ny,2+nz);
b.points(:,3) = b.points(:,3) - 1;
b.points(:,1) = b.points(:,1) - 1;
for n = 1:Nx-2
    b2 = block(Nx-n,Ny,2+nz);
    b2.points(:,3) = b2.points(:,3) - 1 + n*(nz+1);
    b2.points(:,1) = b2.points(:,1) - 1 + n;
    b = mergeBlocks(b,b2,prec);
end

P = b.points;
x = [0 0 0];
n = [-(nz+1) 0 1];
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