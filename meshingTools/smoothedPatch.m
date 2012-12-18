function [x,y] = smoothedPatch(p1,p2,p3,p4,nx,ny,tol)

%p1-p4 should traverse the square in the counterclockwise direction

l1 = lineSegment(p1,p2,nx); 
l2 = lineSegment(p2,p3,ny); 
l3 = lineSegment(p4,p3,nx); 
l4 = lineSegment(p1,p4,ny); 

x = [l1(:,1); l2(:,1); l3(:,1); l4(:,1)];
y = [l1(:,2); l2(:,2); l3(:,2); l4(:,2)];
p = [p1; p2; p3; p4];
xmin = min(p(:,1));
xmax = max(p(:,1));
ymin = min(p(:,2));
ymax = max(p(:,2));
x = xmin + (xmax-xmin)*[0:nx]/nx;
y = ymin + (ymax-ymin)*[0:ny]/ny;
[x,y] = meshgrid(x,y);
% figure(1); clf
% plot(x,y,'.-',x',y','.-r')
% axis equal
%bottom
x(1,:) = l1(:,1)'; 
y(1,:) = l1(:,2)';
%top
x(end,:) = l3(:,1)'; 
y(end,:) = l3(:,2)';
%left
x(:,1) = l4(:,1);
y(:,1) = l4(:,2);
%right
x(:,end) = l2(:,1)';
y(:,end) = l2(:,2)';

% figure(1); hold on
% plot(x,y,'.-',x',y','.-r')
% axis equal
[x,y] = smoothen(x,y,tol);
x = x'; y = y';