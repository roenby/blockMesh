function [F,Ci] = cells2Faces(C)

%[F,Ci] = cells2Faces(C) generates a (nonunique) list of faces F from a 
%list of cells C. C is an Ncells by 8 matrix where the n'th row contains
%8 integers that are indices in a point list of the points comprising the
%corners of the (cubic) cell. Each cell has 6 faces each spanned by 4
%points. Thus, the face list F is a 6*Ncells by 4 matrix and the 4 integers
%in a row are the indices in the point list of the points defining the
%face. By convention, the six faces of a cell are made up of the points 
%with indices:
% I = [ [1 5 8 4]
%       [2 3 7 6]
%       [1 2 6 5]
%       [3 4 8 7]
%       [4 3 2 1]
%       [5 6 7 8] ];

%The second output, Ci, is a list of length 6*Ncells such that the n'th
%element is the index of the cell to which the face belongs. 
%In OpenFOAM speak te face F(n,:) has cell Ci(n) as owner, and therefore 
%the normal defined by the face orientation points out of the cell
%
%Johan Roenby, DHI Water & Environment

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
for n = 1:nc
    c = C(n,:);
    F((1:nf) + (n-1)*nf,:) = c(I);

end
Ci = ceil((1:nf*nc)'/nf);