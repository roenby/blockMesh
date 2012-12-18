function b = unitBlock(nx,ny,nz)

%b = unitBlock(nx,ny,nz) generates a block with nx times ny times nz cells
%going from 0 to 1 in all 3 directions

b = block(nx,ny,nz);
for n = 1:3
    b.points(:,n) = b.points(:,n)-min(b.points(:,n));
    b.points(:,n) = b.points(:,n)/max(b.points(:,n));
end