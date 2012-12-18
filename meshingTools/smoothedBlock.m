function [x,y] = smoothedBlock(l1,l2,l3,l4,tol)

nx = length(l1);
ny = length(l3);
x = [l1(:,1); l2(:,1); l3(:,1); l4(:,1)];
y = [l1(:,2); l2(:,2); l3(:,2); l4(:,2)];
xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);
x = xmin + (xmax-xmin)*[0:nx-1]/(nx-1);
y = ymin + (ymax-ymin)*[0:ny-1]/(ny-1);
[x,y] = meshgrid(x,y);
x = 0*x + 1/2*(xmin + xmax);
y = 0*y + 1/2*(ymin + ymax);
%left side:
x(1,:) = l1(:,1)'; 
y(1,:) = l1(:,2)';
%right side:
x(end,:) = l2(:,1)';
y(end,:) = l2(:,2)';
%bottom
x(:,1) = l3(:,1);
y(:,1) = l3(:,2);
%top
x(:,end) = l4(:,1);
y(:,end) = l4(:,2);

[x,y] = smoothen(x,y,tol);
