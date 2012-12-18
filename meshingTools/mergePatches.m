function b = mergePatches(b,ind,patchName,patchType)

%b = mergePatches(b,ind,patchName,patchType) merges the patches with
%indices ind in the block b and gives the resulting patch the name
%patchName and type patchType
%
%Johan Roenby, DHI Water & Environment

bf1 = b.boundary.startFace(1);
bi = [bf1:b.nFaces];
F = b.faces(bi,:);
O = b.owner(bi);
startFace = b.boundary.startFace;
nFaces = b.boundary.nFaces;


I1 = []; I2 = [];
nn = 0;
for n = 1:length(nFaces)
    fi = startFace(n) - bf1 + [1:nFaces(n)];
    if ismember(n,ind)
        I2 = [I2, fi];
    else
        nn = nn + 1;
        I1 = [I1 fi];
        nFacesNew(nn) = nFaces(n);
        pName{nn} = b.boundary.patchName{n};
        pType{nn} = b.boundary.patchType{n};
    end
end
nFacesNew(nn+1) = length(I2);
startFaceNew = [bf1 + cumsum([0; nFacesNew(1:end-1)'])];
b.faces(bi,:) = F([I1, I2],:);
b.owner(bi) = O([I1, I2]);
b.boundary.nFaces = nFacesNew;
b.boundary.startFace = startFaceNew;

pName{nn+1} = patchName;
pType{nn+1} = patchType;
b.boundary.patchName = pName;
b.boundary.patchType = pType;
