function ind = patchesInPlane(b,x,n,prec)

P = b.points;
bf1 = b.boundary.startFace(1);
bi = [bf1:b.nFaces];
F = b.faces(bi,:);
O = b.owner(bi);
startFace = b.boundary.startFace;
nFaces = b.boundary.nFaces;

[Cf,Sf] = faceCentresAndNormals(P,F);
%quiver3(Cf(:,1),Cf(:,2),Cf(:,3),Sf(:,1),Sf(:,2),Sf(:,3))
absSf = sqrt(sum(Sf.^2,2));
n2 =  Sf./repmat(absSf,1,size(Sf,2)); 

n1 = n(:)/norm(n);
SfEqNormal = abs(n2*n1-1) < prec;
X = repmat(x(:)',size(Cf,1),1);
CfInPlane = abs((Cf-X)*n1) < prec;
isInPlane = SfEqNormal & CfInPlane;
ind = [];
for n = 1:length(startFace)
    fi = startFace(n) - bf1 + [1:nFaces(n)];
    if sum(isInPlane(fi)) == length(fi)
        ind = [ind, n];
    end
end