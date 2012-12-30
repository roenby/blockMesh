function RodsandFluidDomain(caseDir,toolboxDir)

% Generating block meshing of fluid domain around gravity based foundation 
% of the type at the Rodsand II offshore wind farm.
%
% Johan Roenby, DHI Water & Environment

if nargin < 1
    caseDir = pwd;
end
if nargin < 2
    toolboxDir = ['..' filesep '..' filesep 'meshingTools'];
end

%Making case file structure and copying generating m-files to case directory
meshDir = makeCaseDir(caseDir,toolboxDir);
%Copy generating code to case dir
copyGeneratingCode(meshDir,toolboxDir,mfilename('fullpath'));


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

ds = 500; %Basic cell size parameter

zs = cumsum([0, z1, z2, z3, z4, z5]);
z = 0;
for n = 1:length(zs)-1
    nz = round(3*(zs(n+1)-zs(n))/ds);
    z = [z, zs(n) + (zs(n+1)-zs(n))*[1:nz]/nz];
end

prec = 1e-6; %Precision in determining whether two points are identical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defining points

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

r1 = d2/2; r2 = d1/2;
r4 = 20*r1;
r6 = l1/2;
p1 = [d4/2 0];
p2 = [0 d4/2];
p3 = [r2 0];
p4 = [0 r2];
p5 = [r2 r6];
p6 = [r6 r2];

r3 = 2.5*r2;
p7 = [r3 0];
p8 = [r3 r3];
p9 = [0 r3];
p10 = [r4 0];
p11 = [r4 r3];

fac = .2;
p12 = fac*p8+(1-fac)*p9;
p13 = fac*p8+(1-fac)*p7;
p14 = [p10(1) p13(2)];

figure(1); clf
n = 0;
while exist(['p' int2str(n+1)],'var')
    x = eval(['p' int2str(n+1) '(1)']);
    y = eval(['p' int2str(n+1) '(2)']);
    plot(x,y,'.')
    hold on
    text(x,y,['p' int2str(n+1)])
    n = n + 1;
end    
axis equal

%B1
%nPent = round(norm(p13-p8)/2/ds);
%nPent = round(norm(p5-p13)/2/ds);
ny = round(norm(p7-p8)/ds);
nPent = round(ny*norm(p13-p8)/norm(p7-p8)/2);
b1 = smoothedPentagonBlock(p5,p13,p8,p12,p6,nPent,z,prec);

%B2
%n2 = round(norm(p13-p7)/ds);
n2 = ny-2*nPent;
[x,y] = smoothedPatch(p3,p7,p13,p5,2*nPent,n2,prec);
b2 = repeat2DMesh(x,y,z);
b = mergeBlocks(b1,b2,prec);

%B3
[x,y] = smoothedPatch(p4,p6,p12,p9,n2,2*nPent,prec);
b2 = repeat2DMesh(x,y,z);
b = mergeBlocks(b,b2,prec);

%B6
%Inner arc piece
th = lineSegment(0,pi/2,2*n2+2*nPent);
L1 = norm(p1)*[cos(th(:)) sin(th(:))];
%Outer rim
L2a = lineSegment(p3,p5,n2); 
L2b = lineSegment(p5,p6,2*nPent);
L2c = lineSegment(p6,p4,n2);
L2 = [L2a; L2b(2:end,:); L2c(2:end,:)];
n3 = round(2.5*norm(L1(1,:)-L2(1,:))/ds);
[x,y] = rayPatch(L1,L2,n3,1.05);
b2 = repeat2DMesh(x,y,z(z >= zs(2)));
b = mergeBlocks(b,b2,prec);

writePolyMesh(b,meshDir,writePrec,compress)

%Cone
nz = sum(z >= zs(3) & z <= zs(4));
b2 = wedgeBlock(nz-1,2*n2+2*nPent);
R = b2.points(:,1);
TH = b2.points(:,2);
Z = b2.points(:,3);
R = R-min(R); R = R/max(R); R = d2/2 + R*(d4/2-d2/2);
TH = TH-min(TH); TH = TH/max(TH); TH = TH*pi/2;
Z = Z-min(Z); Z = Z/max(Z); Z = zs(3) + (zs(4)-zs(3))*Z;
b2.points(:,1) = R.*cos(TH); 
b2.points(:,2) = R.*sin(TH); 
b2.points(:,3) = Z; 
b = mergeBlocks(b,b2,prec);

%Middle cylinder
L0 = r1*[cos(th(:)) sin(th(:))];
%n4 = round(norm(L0(1,:)-L1(1,:))/ds);
[x,y] = rayPatch(L0,L1,nz);
b2 = repeat2DMesh(x,y,z(z >= zs(2) & z <= zs(3)));
b = mergeBlocks(b,b2,prec);

%Upper cylinder
[x,y] = rayPatch(L0,L1,ceil(.8*nz),1.02); %Refining in radial direction due to crashing solution
b2 = repeat2DMesh(x,y,z(z >= zs(5) & z <= zs(6)));
b = mergeBlocks(b,b2,prec);

%Copying to second quadrant
b2 = rotBlock(b,[0 0 pi/2]);
b = mergeBlocks(b,b2,prec);

%B4
n = round(norm(p7-p10)/ds);
[x,y] = smoothedPatch(p7,p10,p14,p13,n,n2,prec);
b2 = repeat2DMesh(x,y,z);

%B5
[x,y] = smoothedPatch(p13,p14,p11,p8,n,2*nPent,prec);
b3 = repeat2DMesh(x,y,z);
b2 = mergeBlocks(b2,b3,prec);

b = mergeBlocks(b,b2,prec);

x = b2.points(:,1);
x = x - min(x) - max(x);
b2.points(:,1) = x;
b = mergeBlocks(b,b2,prec);

%Making other half of domain
b2 = rotBlock(b,[0 0 pi]);
b = mergeBlocks(b,b2,prec);

%Merging patches
ind = patchesInPlane(b,[0 -p9(2) 0],[0 -1 0],prec);
b = mergePatches(b,ind,'front','patch');

ind = patchesInPlane(b,[0 p9(2) 0],[0 1 0],prec);
b = mergePatches(b,ind,'back','patch');

ind = patchesInPlane(b,[0 0 zs(end)],[0 0 1],prec);
b = mergePatches(b,ind,'top','patch');

ind = patchesInPlane(b,[0 0 zs(1)],[0 0 -1],prec);
b = mergePatches(b,ind,'bottom','patch');

ind = patchesInPlane(b,[p10(1) 0 0],[1 0 0],prec);
b = mergePatches(b,ind,'outlet','patch');

ind = patchesInPlane(b,[-p10(1) 0 0],[-1 0 0],prec);
b = mergePatches(b,ind,'inlet','patch');

patchInd = strfind(b.boundary.patchName,'patch');
ind = zeros(length(patchInd),1);
for n = 1:length(patchInd)
    if ~isempty(patchInd{n})
        ind(n) = 1;
    end
end
ind = find(ind);
b = mergePatches(b,ind,'structure','patch');

%Move z-level to water surface
%waterLevelInd = find(abs(z-waterLevel) == min(abs(z-waterLevel)));
%b.points(:,3) = b.points(:,3) - z(waterLevelInd);

%Rescale to meters
fac = 1000;
b.points = b.points/fac;

writePolyMesh(b,meshDir,writePrec,compress)

%Write z-levels to file (to be used when using water level with setFields)
z = z/fac;
fid = fopen([meshDir filesep 'zLevels'],'w');
nsd = ceil(log10(max(abs(z)))) + writePrec; %number of significant digits
formatString = ['%0.' int2str(nsd) 'g\n'];
fprintf(fid, formatString, z);
fclose(fid);
