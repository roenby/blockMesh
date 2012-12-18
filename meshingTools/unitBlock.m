function b = unitBlock(nx,ny,nz)

b = block(nx,ny,nz);
for n = 1:3
    b.points(:,n) = b.points(:,n)-min(b.points(:,n));
    b.points(:,n) = b.points(:,n)/max(b.points(:,n));
end