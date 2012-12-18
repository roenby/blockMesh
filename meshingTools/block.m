function b = block(nx,ny,nz)

% b = block(nx,ny,nz) generates a block mesh with grid points [0:nx] along
% the x-axis, [0:ny] along the y-axis, and [0:nz] along the z-axis. 
% The structure returned contains the necessary data to generate the
% files points, faces, owner, neighbour, and boundary that define a
% polyhedral mesh in OpenFOAM:
%
% b.points = List of grid points, (nPoints x 3)-matrix
% b.faces = List of faces, (nFaces x 4)-matrix with indices in b.points
% b.owner = Array of cells to which the faces "belong", length = nFaces.
% b.neighbour = Array of other cell that shares the face with the owner cell, length = nInternalFaces
%b.boundary.startFace = Array containing the start indices in the face list for each of the 6 boundary patches
%b.boundary.nFaces = Array containing the number of faces for each of the 6 boundary patches.
%b.boundary.patchName = {'left' 'right' 'front' 'back' 'bottom' 'top'};
%b.boundary.patchType = {'patch' 'patch' 'patch' 'patch' 'patch' 'patch'};
%
% Any connectivity and topology preserving transformation may be performed 
% on the points without any need for altering the remainin data.
%
% Johan Roenby, DHI Water & Environment

if nargin < 1
    nx = 5;
    ny = 7;
    nz = 6;
end

%Changing to new convention that input to this function is the number of 
%bins/cells in each direction.
nx = nx + 1; ny = ny + 1; nz = nz + 1;
%From here on nx, ny and nz denote the number of grid points rather than bins.

nPoints = nx*ny*nz;
nCells = (nx-1)*(ny-1)*(nz-1);
nBoundaryFaces = 2*((nx-1)*(ny-1) + (nx-1)*(nz-1) + (ny-1)*(nz-1));
nUniqueFaces = nx*(ny-1)*(nz-1) + ny*(nz-1)*(nx-1) + nz*(nx-1)*(ny-1);
nInternalFaces = nUniqueFaces - nBoundaryFaces;

b.nPoints = nPoints;
b.nCells = nCells;
b.nInternalFaces = nInternalFaces;
b.nFaces = nUniqueFaces;

%Generating cell point index list
C = cellIndexList(nx,ny,nz);

%Old plotting stuff
% P = [Xm(:) Ym(:) Zm(:)]
% figure(1); clf
% m = 1;
% for n = 1:8
%     pmn = P(C(m,n),:);
%     plot3(pmn(1),pmn(2),pmn(3),'.')
%     hold on
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
% end

%The six face of a cell are made up of the points with indices:
[F,Ci] = cells2Faces(C);

%Old plotting stuff
% for n = 1:size(F,1)
%     fp = P(F(n,:),:);
%     cp = P(C(Ci(n),:),:);
%     Co = mean(cp,1);
%     Cf = mean(fp,1);
%     Sf = cross(fp(2,:)-fp(1,:),fp(3,:)-fp(2,:));
%     if dot((Cf-Co),Sf) < 0
%         disp(['Face ' int2str(n) ' has wrong orientation'])
%     end
% end

%Find pairs of faces that are the same face seen from the owner and the
%neighbour cells
[~,I1] = unique(sort(F,2),'rows','first');
[~,I2] = unique(sort(F,2),'rows','last');

%I1(I1~=I2) are the indices in the face list of the first apperance of all 
%face appearing twice in the face list, i.e. internal faces. 
%Ci(I1(I1~=I2)) is therefore the corresponding owner cell list.
%I2(I1~=I2) are the indices in the face list of the second apperance of all
%faces appearing twice in the face list, i.e. internal faces. 
%Ci(I2(I1~=I2)) is therefore the corresponding neighbour cell list.
%I1(I1==I2) are the indices in the face list of the faces only appearing 
%once, i.e. boundary faces.
%Ci(I1(I1==I2)) is therefore the corresponding owner cell list for boundary
%faces.

%Removing second apperance of internal faces in face list and moving 
%boundary faces to the back of the face list

F = F([I1(I1~=I2); I1(I1==I2)],:);
owner = Ci([I1(I1~=I2); I1(I1==I2)]);
neighbour = Ci(I2(I1 ~= I2));

%Sorting internal part of face list by owner and neighbour cell 
[~,I] = sortrows([owner(1:nInternalFaces) neighbour]);
F(1:nInternalFaces,:) = F(I,:);
owner(1:nInternalFaces) = owner(I);
neighbour = neighbour(I);

%Old plotting stuff
%ind = F(nInternalFaces+1:end,:);
%plot3(Xm(ind),Ym(ind),Zm(ind),'.')
%Seems to get the boundary points correct

%Sorting boundary part of face list by side and owner
nFaces = [  (ny-1)*(nz-1); (ny-1)*(nz-1); 
            (nx-1)*(nz-1); (nx-1)*(nz-1);
            (nx-1)*(ny-1); (nx-1)*(ny-1);];

startFace = [nInternalFaces + 1 + cumsum([0; nFaces(1:end-1)])];

clear I
[I{1},I{2},I{3}] = ind2sub([nx ny nz],F(nInternalFaces+1:end,:));
lim = [1 nx; 1 ny; 1 nz];
ind = [];%zeros(sum(nFaces),1);
for i = 1:length(I)
    for j = 1:size(lim,2)
        ind_ij = sort(find(sum(I{i},2) == 4*lim(i,j)));
        ind = [ind; ind_ij];
    end
end
ind = ind + nInternalFaces;
F(nInternalFaces+1:end,:) = F(ind,:);
owner(nInternalFaces+1:end) = owner(ind);

%Generating grid points
[Xm,Ym,Zm] = ind2sub([nx ny nz],1:nx*ny*nz);

%Putting data into structure
b.points = [Xm(:)-1 Ym(:)-1 Zm(:)-1]; %Changed so it goes from 0 to nx, ny and nz instead of from 1
b.faces = F;
b.owner = owner;
b.neighbour = neighbour;
b.boundary.nFaces = nFaces;
b.boundary.startFace = startFace;

b.boundary.patchName = {'left' 'right' 'front' 'back' 'bottom' 'top'};
b.boundary.patchType = {'patch' 'patch' 'patch' 'patch' 'patch' 'patch'};