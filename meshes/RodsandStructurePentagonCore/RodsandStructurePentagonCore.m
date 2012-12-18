function RodsandStructurePentagonCore(caseDir)

%Generating structural block mesh for gravity based foundation of the type 
%installed at the Rodsand 2 offshore farm.
%
%The core of the structure is filled with a pentagon.
%
%Johan Roenby, DHI Water & Environment

if nargin < 1
    caseDir = pwd;
end

%Making case file structure and copying generating m-files to case directory
toolboxDir = ['..' filesep '..' filesep 'meshingTools'];
copyfile(toolboxDir,[caseDir filesep 'private'])
meshDir = makeCaseDir(caseDir);

compress = 1; %1 to compress output files 
writePrec = 12; %Write precision in output files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defining dimensions

z1 = 3450; %Foot height
z2 = 2650; %Middle cylinder height
z3 = 3800; %Cone height
z4 = 2100; %Top cylinder height
z5 = z3/2; %Tower section height
l1 = 4650; %Short side length
d1 = 17000; %Foot diameter at short side
d2 = 4300; %Middle cylinder diameter
d4 = 9000; %Top cylinder diameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defining grid control parameters

ds = 250; %Basic cell size parameter

zs = cumsum([0, z1, z2, z3, z4, z5]);
z = 0;
for n = 1:length(zs)-1
    nz = round(1*(zs(n+1)-zs(n))/ds);
    z = [z, zs(n) + (zs(n+1)-zs(n))*[1:nz]/nz];
end

%Number of azimuthal bins - MUST BE DIVISIBLE BY 4 AND 10!!!
nth = 10*8;

%Inner core parameters
nr2 = 8; %Number of radial bins of non-pentagon parts
r2 = d2/2; %Core radius
r1 = .75*r2; %Pentagon radius

prec = 1e-6; %Precision in determining whether two points are identical

%Rescaling factor for mesh
fac = 1e-3; %mm -> m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defining points for fluid domain

%Layout of blocks as seen from above:
%
%p9 ______p12____________p8______________________________p11
%  |      |              |                               |
%  |  B3  |              |                               |
%p4|______|p6     B1     |               B5              |
%  |       \             |                               |
%  |        \            |                               |
%  |__   B6  \p5_________|p13____________________________|p14
%p2    \      |    B2    |               B4              |
%       |_____|__________|_______________________________|
%       p1    p3         p7                             p10

p3 = [d1/2 0];
p4 = [0 d1/2];
p5 = [d1/2 l1/2];
p6 = [l1/2 d1/2];

%Foot
th = lineSegment(0,pi/2,nth/4);
L1 = d2/2*[cos(th(:)) sin(th(:))];
nth1 = round(nth/4*abs(p3-p5)/(2*abs(p3-p5) + abs(p5-p6)));
nth2 = nth/4 - 2*nth1;
L2a = lineSegment(p3,p5,nth1); 
L2b = lineSegment(p5,p6,nth2);
L2c = lineSegment(p6,p4,nth1);
L2 = [L2a; L2b(2:end,:); L2c(2:end,:)];
n3 = round(norm(L1(1,:)-L2(1,:))/ds);
[x,y] = rayPatch(L1,L2,n3,1.1);
b = repeat2DMesh(x,y,z(z <= zs(2)));
b2 = b;
for n = 1:3
    b2 = rotBlock(b2,[0 0 pi/2]);
    b = mergeBlocks(b,b2,prec);
end

%Inner core
b2 = pentagonInCircle(nth/10,nr2,r1,r2,z,prec);
b = mergeBlocks(b,b2,prec);

%Cone
nz = sum(z >= zs(3) & z <= zs(4));
b2 = wedgeBlock(nz,nth/2+1,prec);
b2 = rotBlock(b2,[0 pi 0]);
b2.points(:,1) = b2.points(:,1)-min(b2.points(:,1));
b2.points(:,2) = b2.points(:,2)-min(b2.points(:,2));
b2.points(:,3) = b2.points(:,3)-min(b2.points(:,3));
R = b2.points(:,1);
TH = b2.points(:,2);
Z = b2.points(:,3);
R = R-min(R); R = R/max(R); R = d2/2 + R*(d4/2-d2/2);
TH = TH-min(TH); TH = TH/max(TH); TH = TH*pi;
Z = Z-min(Z); Z = Z/max(Z); Z = zs(3) + (zs(4)-zs(3))*Z;
b2.points(:,1) = R.*cos(TH); 
b2.points(:,2) = R.*sin(TH); 
b2.points(:,3) = Z; 
b3 = mergeBlocks(b2,rotBlock(b2,[0 0 pi]),prec);
b = mergeBlocks(b,b3,prec);

%Upper cylinder rim
nz2 = sum(z >= zs(4) & z <= zs(5));
b2 = unitBlock(nz,nth/2+1,nz2);
R = b2.points(:,1);
TH = b2.points(:,2);
Z = b2.points(:,3);
R = R-min(R); R = R/max(R); R = d2/2 + R*(d4/2-d2/2);
TH = TH-min(TH); TH = TH/max(TH); TH = TH*pi;
Z = Z-min(Z); Z = Z/max(Z); Z = zs(4) + (zs(5)-zs(4))*Z;
b2.points(:,1) = R.*cos(TH); 
b2.points(:,2) = R.*sin(TH); 
b2.points(:,3) = Z; 
b3 = mergeBlocks(b2,rotBlock(b2,[0 0 pi]),prec);
b = mergeBlocks(b,b3,prec);

%Merging patches
ind = patchesInPlane(b,[0 0 zs(end)],[0 0 1],prec);
b = mergePatches(b,ind,'top','patch');

ind = patchesInPlane(b,[0 0 0],[0 0 -1],prec);
b = mergePatches(b,ind,'bottom','patch');

patchInd = strfind(b.boundary.patchName,'patch');
ind = zeros(length(patchInd),1);
for n = 1:length(patchInd)
    if ~isempty(patchInd{n})
        ind(n) = 1;
    end
end
ind = find(ind);
b = mergePatches(b,ind,'theRest','patch');

%Rescaling and printing mesh
b.points = fac*b.points;

writePolyMesh(b,meshDir,writePrec,compress)

%plotMesh(b)