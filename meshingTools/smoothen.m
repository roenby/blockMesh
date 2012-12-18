function [x,y] = smoothen(x,y,tol)

%[x,y] = smoothen(x,y,tol) smoothens the positions of all non-boundary 
%points of the mesh (x,y). In the resulting mesh the coordinates of an 
%internal point will be the average of the nearest neighbour points to
%within the specified tolerance tol.
%
%Johan Roenby, DHI Water & Environment

err = 2*tol;
while err > tol
    xold = x;
    yold = y;
    s1 = 1;
    x(2:end-1,2:end-1) = ( ...
        s1*x(1:end-2,2:end-1) + (2-s1)*x(3:end,2:end-1) + ...
        x(2:end-1,1:end-2) + x(2:end-1,3:end)...
        )/4;
    y(2:end-1,2:end-1) = ( ...
        y(1:end-2,2:end-1) + y(3:end,2:end-1) + ...
        y(2:end-1,1:end-2) + y(2:end-1,3:end)...
        )/4;
    err = max(hypot(x(:)-xold(:),y(:)-yold(:)));
%     plot(x,y,'r.')
%     axis equal
%     drawnow
%     hold off
end