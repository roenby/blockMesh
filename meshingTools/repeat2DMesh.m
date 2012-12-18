function b = repeat2DMesh(X,Y,z)

[nx ny] = size(X);
nz = length(z);
X3D = zeros(nx,ny,nz);
Y3D = zeros(nx,ny,nz);
Z3D = zeros(nx,ny,nz);
for n = 1:nz;
    X3D(:,:,n) = X;
    Y3D(:,:,n) = Y;
    Z3D(:,:,n) = z(n)*ones(size(X));
end
b = block(nx,ny,nz);
b.points = [X3D(:) Y3D(:) Z3D(:)];
b.boundary.patchName = {'left' 'right' 'front' 'back' 'bottom' 'top'};
b.boundary.patchType = {'patch' 'patch' 'patch' 'patch' 'patch' 'patch'};
