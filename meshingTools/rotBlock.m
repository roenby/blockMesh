function b = rotBlock(b,v)

P = b.points;
nComp = size(P,2);
th = norm(v);
e1 = repmat(v(:)'/th,size(P,1),1);
e2 = cross(e1,cross(P,e1)); %Component of P perpendicular to v
abse2 = sqrt(sum(e2.^2,2));

%Sorting out points that are parallel to v;
ind1 = find(abse2 ~= 0);
ind2 = find(abse2 == 0);
Pnew = P(ind1,:);
e1 = e1(ind1,:);
e2 = e2(ind1,:);
abse2 = repmat(abse2(ind1),1,nComp);
e2 = e2./abse2;
e3 = cross(e1,e2);

Ppara = dot(Pnew,e1,2);
Pperp = dot(Pnew,e2,2);
e2new = cos(th).*e2 + sin(th).*e3;
P(ind1,:) = repmat(Ppara,1,nComp).*e1 + repmat(Pperp,1,nComp).*e2new;
b.points = P;