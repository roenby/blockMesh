function C = cellIndexList(nx,ny,nz)

%The n'th row in C = cellIndexList(nx,ny,nz) contains the indices of the
%points in the grid point list belonging to the n'th cell of a block with
%(nx x ny x nz) grid points.

I = [   [0 0 0]
        [1 0 0]
        [1 1 0]
        [0 1 0]
        [0 0 1]
        [1 0 1]
        [1 1 1]
        [0 1 1] ];

nCells = (nx-1)*(ny-1)*(nz-1);
nPointsPerCell = size(I,1);
C = zeros(nCells,nPointsPerCell);
for i = 1:nPointsPerCell
    [I2,I1,I3] = meshgrid((1:ny-1)+I(i,2),(1:nx-1)+I(i,1),(1:nz-1)+I(i,3));
    ind = sub2ind([nx ny nz],I1,I2,I3);
    C(:,i) = ind(:);
end
