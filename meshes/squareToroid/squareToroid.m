function squareToroid(caseDir,toolboxDir)

% Generating toroidal (doughnut shaped) block mesh with a square cross
% section.
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


compress = 0; %1 to compress output files 
writePrec = 12; %Write precision in output files
prec = 1e-6; %Precision in determining whether two points are identical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defining dimensions
l = 1; %cross section side length
R = 2.5; %Major radius

nl = 5; %Number of bins along sides
nth = 40; %Number of azimuthal bins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Making half a torus
b = unitBlock(nl,nth/2,nl);
rr = R + l*(b.points(:,1)-.5);
th = pi*b.points(:,2);
Z = l*(b.points(:,3)-.5);
X = rr.*cos(th);
Y = rr.*sin(th);
b.points(:,1) = X;
b.points(:,2) = Y;
b.points(:,3) = Z;

%Copy-pasting to complete torus
b2 = b;
b2 = rotBlock(b2,[0 0 pi]);
b = mergeBlocks(b,b2,prec);

%Mergin patches
ind = patchesInPlane(b,[0 0 l/2],[0 0 1],prec);
b = mergePatches(b,ind,'top','patch');
ind = patchesInPlane(b,[0 0 -l/2],[0 0 -1],prec);
b = mergePatches(b,ind,'bottom','patch');

%plotMesh(b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Writing to OpenFOAM polyMesh files
writePolyMesh(b,meshDir,writePrec,compress)
