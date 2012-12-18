function [F,Ci] = cells2Faces(C)

%The six faces of a cell are made up of the points with indices:
% I = [ [1 5 8 4]
%       [2 3 7 6]
%       [1 2 6 5]
%       [3 4 8 7]
%       [4 3 2 1]
%       [5 6 7 8] ];

%Arranged such that F(n,:) has Ci(n) as owner, i.e. the normal defined by
%the face orientation points out of the cell
I = [ [4 3 2 1]
      [5 6 7 8]
      [1 2 6 5]
      [3 4 8 7]
      [1 5 8 4]
      [2 3 7 6] ];

  
%Calculating faces from cells
nf = size(I,1); %Number of faces per cell
nc = size(C,1); %Number of cells
np = size(I,2); %Number of points per face
F = zeros(nf*nc,np); %Face list with internal faces included twice
%Ci = zeros(nc*nf,1); %Cell index list
%fn = 0;
for n = 1:nc
%    ind = (1:nf) + (n-1)*nf;
%    c = C(n,:);
%    F(ind,:) = c(I);
%    ind = (1:nf) + (n-1)*nf;
    c = C(n,:);
    F((1:nf) + (n-1)*nf,:) = c(I);

%    Ci(ind) = n*ones(nf,1);
%    for nf = 1:size(I,1)
%        fn = fn + 1;
%        F(fn,:) = c(I(nf,:));
%        Ci(fn) = nb;
%    end
end
Ci = ceil((1:nf*nc)'/nf);