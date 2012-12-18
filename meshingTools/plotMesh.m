function plotMesh(b)

P = b.points;
plot3(P(:,1),P(:,2),P(:,3),'.k','markersize',2)
hold on
axis equal
cols = 'rgbkcmyrgbkcmyrgbkcmyrgbkcmy';
for n = 1:size(b.faces,1)
    if n == length(b.neighbour)+1
%        clf
        view(3)
        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
    P = b.points(b.faces(n,:),:);
    Co = calcCellCentre(b.owner(n),b);
    Cf = mean(P,1);
    Sf = cross(P(2,:)-P(1,:),P(3,:)-P(2,:));
    if dot((Cf-Co),Sf) < 0
        disp(['Face ' int2str(n) ' has wrong orientation'])
    end
    col = cols(find(b.boundary.startFace > n,1));
    h = patch(P(:,1),P(:,2),P(:,3),col);
    hold on
    plot3(Co(1),Co(2),Co(3),'.g')
    if 0 %n <= length(b.neighbour)
        Cn = calcCellCentre(b.neighbour(n),b);
        C = [Co(:) Cn(:)]';
        plot3(C(:,1),C(:,2),C(:,3),'k')
        plot3(Cn(1),Cn(2),Cn(3),'oc')
        quiver3(Cf(1),Cf(2),Cf(3),Sf(1),Sf(2),Sf(3))
        plot3(P(:,1),P(:,2),P(:,3),'*g')
        xlabel('x')
        ylabel('y')
        zlabel('z')
    end
    set(h,'facealpha',.2)
    drawnow
    set(h,'facealpha',.5)
end

function C = calcCellCentre(n,b)

ownerFaceInd = find(b.owner == n);
neighbourFaceInd = find(b.neighbour == n);
nFaces = length(ownerFaceInd) + length(neighbourFaceInd);
if nFaces ~= 6
    disp(['Cell ' int2str(n) ' has ' int2str(nFaces) ' faces!!!!!!!'])
    disp('Faces owned:')
    disp(num2str(ownerFaceInd))
    disp('Faces neighboured:')
    disp(num2str(neighbourFaceInd))
end
edgeInds = unique(b.faces([ownerFaceInd(:); neighbourFaceInd(:)],:));
C = sum(b.points(edgeInds,:),1)/numel(edgeInds);