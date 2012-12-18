function [x,y] = smoothedPatch(p1,p2,p3,p4,nx,ny,tol)

%[x,y] = smoothedPatch(p1,p2,p3,p4,nx,ny,tol) generates a 2d mesh between
%the 4 points p1-p4, which should traverse the square in the 
%counterclockwise direction. There will be nx cells or bins on hte line 
%from p1 to p2 and from p3 to p4. There will be ny points from p2 to p3 and
%from p4 to p1. The inner points are smoothely distributed.
%
%Johan Roenby, DHI Water & Environment

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