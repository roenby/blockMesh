function b = mergeBlocks(b1,b2,prec)

nif1 = b1.nInternalFaces;
nif2 = b2.nInternalFaces;
nf1 = b1.nFaces;
nf2 = b2.nFaces;
np1 = b1.nPoints;
np2 = b2.nPoints;
nc1 = b1.nCells;
nc2 = b2.nCells;
P1 = b1.points;
P2 = b2.points;
F1 = b1.faces;
F2 = b2.faces;
O1 = b1.owner;
O2 = b2.owner;
N1 = b1.neighbour;
N2 = b2.neighbour;

%FINDING IDENTICAL BOUNDARY POINTS IN b1 AND b2

bf1 = F1(nif1+1:end,:); %Boundary faces in b1
Ibp1 = unique(bf1(:)); %Indices of boundary points in b1.points
bp1 = P1(Ibp1,:); %Boundary points from b1.points
bf2 = F2(nif2+1:end,:); %Boundary faces in b2
Ibp2 = unique(bf2(:));%Indices of boundary points in b2.points
bp2 = P2(Ibp2,:); %Boundary points from b2.points

%I1 and I2 are the indices in the boundary part of b1.points and b2.points 
%that are the same to within the precision prec.
[~,I1,I2] = intersect(round(bp1/prec)*prec,round(bp2/prec)*prec,'rows');
%The corresponding point indices in b1.points and b2.points are:
Ip1 = Ibp1(I1); Ip2 = Ibp2(I2); 

% figure(1); clf
% plot3(bp1(:,1),bp1(:,2),bp1(:,3),'.')
% hold on 
% plot3(bp2(:,1),bp2(:,2),bp2(:,3),'or')
% figure(1); clf
% plot3(b1.points(Ip1,1),b1.points(Ip1,2),b1.points(Ip1,3),'oc')
% hold on 
% plot3(b2.points(Ip2,1),b2.points(Ip2,2),b2.points(Ip2,3),'.y')

%Indices of point in b2.points that are not in b1:
notIp2 = setdiff(1:np2,Ip2);
%Making a point list for the merged blocks containing all the points
%from b1 and all the points in b2 that are not in b1
b.points = [P1; P2(notIp2,:)];
%The points b2.points(Ip2,:) now have indices Ip1 in b.points while the 
%points b2.points(notIp2,:) have indices:
notIp1 = (1:length(notIp2)) + np1;
%Mapping of boundary points indices in b2 to indices in b
%bpMap = sortrows([Ip2(:), Ip1(:); notIp2(:), notIp1(:)]);
bpMap = sortrows([ [Ip2(:); notIp2(:)],  [Ip1(:);  notIp1(:)] ]);
%Applying the map so that pointers to points in b2 that are also in b1 are
%changed to point to the b1 version of the point and pointers to unique
%boundary points in b2 are renumbered to correspond to their indices in
%the new b.points point array
b2Faces = interp1(bpMap(:,1),bpMap(:,2),F2);
%Now we can compare the boundary faces of b1 and b2 via the indices
b1bFaces = F1(nif1+1:nf1,:);
b2bFaces = b2Faces(nif2+1:nf2,:);
[~,I1,I2] = intersect(sort(b1bFaces,2),sort(b2bFaces,2),'rows');
%I1 and I2 are index arrays of the same length such that the face
%b2bFaces(I2(n),:) is identical to the face b1bFaces(I1(n),:)
%The newly created internal faces are now
newInternalFaces = b1bFaces(I1,:);
%We now remove boundary faces that have become internal faces in the merged 
%mesh from the boundary face lists:
notI1 = setdiff(1:size(b1bFaces,1),I1); %Indices of boundary faces in b1bFaces that are not in b2bFaces 
b1bFaces = b1bFaces(notI1,:);
notI2 = setdiff(1:size(b2bFaces,1),I2); %Indices of boundary faces in b2bFaces that are not in b1bFaces 
b2bFaces = b2bFaces(notI2,:);

%The merged face list is now
b.faces = [ F1(1:nif1,:);        %Internal faces from b1
            newInternalFaces;    %New internal faces 
            b2Faces(1:nif2,:);   %Internal faces from b2
            b1bFaces;             %boundary faces from b1
            b2bFaces            %boundary faces from b2
           ];

%In the new merged mesh we make the convention that the cells from b2 are
%appended to the cells from b1. Thus the cell numbers of the b1 cells
%remain the same while b1.nCells should be added to the cell numbers of the
%b2 cells.

%The internal faces of b1 therefore keep their original owner and neighbour

%The new internal faces also keep their owner
b1bFaceOwner = O1(nif1+1:nf1);
newInternalFaceOwner = b1bFaceOwner(I1);
%They did not have a neighbour in the b1 mesh but in the merged mesh their
%neighbours are the owners of the corresponding boundary face in b2:
b2bFaceOwner = O2(nif2+1:nf2);
newInternalFaceNeighbour = b2bFaceOwner(I2) + nc1;

%The internal faces of b2 keep their original owner and neighbour cells
%except that b1.nCells should be added to these cell numbers
b2iFaceOwner = O2(1:nif2) + nc1;
b2iFaceNeighbour = N2(1:nif2) + nc1;

%The b1 and b2 boundary face owner lists should be cleaned so they no 
%longer contain the entries corresponding to the new internal faces
b1bFaceOwner = O1(nif1+1:nf1); %Original b1 boundary owner list
b1bFaceOwner = b1bFaceOwner(notI1); %Sorting out new internal face entries
b2bFaceOwner = O2(nif2+1:nf2); %Original b1 boundary owner list
b2bFaceOwner = b2bFaceOwner(notI2) + nc1; %Sorting out new internal face entries

b.owner = [ O1(1:nif1);
            newInternalFaceOwner;
            b2iFaceOwner;
            b1bFaceOwner;
            b2bFaceOwner
           ];
            
b.neighbour = [ N1(1:nif1); newInternalFaceNeighbour; b2iFaceNeighbour];

%plot3(b.points(b2bFaces,1),b.points(b2bFaces,2),b.points(b2bFaces,3),'.r')

%Removing duplicate faces from b2bFaces and making their owners into the
%neighbour of the corresponding face in b1bFaces


% clf
% p = b.points(b1bFaces(notI1,:),:);
% plot3(p(:,1),p(:,2),p(:,3),'k.')
% hold on
% p = b2.points;
% plot3(p(:,1),p(:,2),p(:,3),'or')

%Making boundary info for merged block
nn = 0;
for n = 1:length(b1.boundary.nFaces)
    nFaces = b1.boundary.nFaces(n);
    startFace = b1.boundary.startFace(n);
    ind = startFace + [0:nFaces-1] - nif1;
    nNewIntFaces = sum(I1>=ind(1) & I1<=ind(end));%More efficient than setdiff
    nFaces = nFaces - nNewIntFaces;%
    if nFaces %numel(setdiff(ind,I1))
        nn = nn + 1;
        b.boundary.nFaces(nn) = nFaces; %numel(setdiff(ind,I1));
        b.boundary.patchType{nn} = b1.boundary.patchType{n};
%        b.boundary.patchName{nn} = b1.boundary.patchName{n};
    end
end
for n = 1:length(b2.boundary.nFaces)
    nFaces = b2.boundary.nFaces(n);
    startFace = b2.boundary.startFace(n);
    ind = startFace + [0:nFaces-1] - nif2;
    nNewIntFaces = sum(I2>=ind(1) & I2<=ind(end));%More efficient than setdiff
    nFaces = nFaces - nNewIntFaces;%
    if nFaces %numel(setdiff(ind,I2))
        nn = nn + 1;
        b.boundary.nFaces(nn) = nFaces; %numel(setdiff(ind,I2));
        b.boundary.patchType{nn} = b2.boundary.patchType{n};
%        b.boundary.patchName{nn} = b2.boundary.patchName{n};
    end
end
b.nInternalFaces = nif1 + nif2 + length(I1);
b.boundary.startFace = [b.nInternalFaces + 1 + cumsum([0; b.boundary.nFaces(1:end-1)'])];
b.nPoints = size(b.points,1);
b.nCells = nc1 + nc2;
b.nFaces = nf1 + nf2 - length(I1);

%Ensuring that no two patches has the same name by adding numbers to
%patches with the same name
for n = 1:length(b.boundary.nFaces)
%    ind = find(strcmp(b.boundary.patchName,b.boundary.patchName{n}));
%    ind = ind(ind >= n);
%    if length(ind) > 1
%        for m = 1:length(ind)
%            b.boundary.patchName{ind(m)} = [b.boundary.patchName{ind(m)} int2str(m-1)];
%        end
%    end
     b.boundary.patchName{n} = ['patch' int2str(n)];
end
%Sorting boundary face patches of merged block by ascending owner cell number
% [~,I] = sort(b.owner(b.boundary.startFace));
% bFaces = []; owner = [];
% for n = 1:length(I)
%     startFace = b.boundary.startFace(I(n));
%     nFaces(n) = b.boundary.nFaces(I(n));
%     endFace = startFace + nFaces(n)-1;
%     bFaces = [bFaces; b.faces(startFace:endFace,:)];
%     owner = [owner; b.owner(startFace:endFace)];
% end
% b.boundary.nFaces = nFaces;
% b.boundary.startFace = [b.nInternalFaces + 1 + cumsum([0; b.boundary.nFaces(1:end-1)'])];
% b.faces(b.boundary.startFace(1):end,:) = bFaces;
% b.owner(b.boundary.startFace(1):end) = owner;