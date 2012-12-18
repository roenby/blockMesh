function writePolyMesh(b,polyMeshDir,prec,compress)

%Write block mesh structure to polyMesh files: points, faces, owner,
%neighbour, and boundary.
%
%Note: Unfortunately this becomes rather slow for large meshes.
%
%Johan Roenby, DHI Water & Environment

if nargin < 2
    polyMeshDir = [pwd '\constant\polyMesh'];
end
if nargin < 3
    prec = 12;
end
if nargin < 4
    compress = 0;
end

%nl = char(10); %New line character in Unix ([char(10),char(13)] in windows and char(13) on mac)
footer =  '// ************************************************************************* //';

%chunkSize = 1e5;

fid = fopen('points','w');
writeHeader(fid,'vectorField','points');
fprintf(fid, '%u\n', b.nPoints);
fprintf(fid, '%s\n', '(');
nsd = ceil(log10(max(abs(b.points(:))))) + prec; %number of significant digits
P = round(b.points*10^prec)/10^prec;
formatString = ['(%0.' int2str(nsd) 'g %0.' int2str(nsd) 'g %0.' int2str(nsd) 'g)\n'];
%n1 = 1;
%while n1 <= size(P,1)
%    n2 = min(n1 + chunkSize,size(P,1));
%    fprintf(fid, formatString, P(n1:n2,:)');
%    n1 = n2 + 1;
%end
fprintf(fid, formatString, P');
fprintf(fid, '%s\n', ')');
fprintf(fid, '\n\n%s\n', footer);
fclose(fid);

delete([polyMeshDir filesep 'points*'])
if compress
    gzip('points',polyMeshDir)
    delete('points')
else
    movefile('points',polyMeshDir)
end

%Write faces file
fid = fopen('faces','w');
writeHeader(fid,'faceList','faces');
fprintf(fid, '%u\n', b.nFaces);
fprintf(fid, '%s\n', '(');
nVert = sum(b.faces ~= -1,2);
F = [nVert, b.faces-1];
%n1 = 1;
%while n1 <= size(F,1)
%    n2 = min(n1 + chunkSize,size(F,1));
%    fprintf(fid, formatString, F(n1:n2,:)');
%    n1 = n2 + 1;
%end
fprintf(fid, '%u(%u %u %u %u)\n', F');
fprintf(fid, '%s\n', ')');
fprintf(fid, '\n%s\n', footer);
fclose(fid);

delete([polyMeshDir filesep 'faces*'])
if compress
    gzip('faces',polyMeshDir)
    delete('faces')
else
    movefile('faces',polyMeshDir)
end

%Write owner file
fid = fopen('owner','w');
writeHeader(fid,'labelList','owner',b);
fprintf(fid, '%u\n', b.nFaces);
fprintf(fid, '%s\n', '(');
fprintf(fid, '%u\n', b.owner-1);
fprintf(fid, '%s\n', ')');
fprintf(fid, '\n%s\n', footer);
fclose(fid);

delete([polyMeshDir filesep 'owner*'])
if compress
    gzip('owner',polyMeshDir)
    delete('owner')
else
    movefile('owner',polyMeshDir)
end

%Write neighbour file
fid = fopen('neighbour','w');
writeHeader(fid,'labelList','neighbour',b);
fprintf(fid, '%u\n', b.nInternalFaces);
fprintf(fid, '%s\n', '(');
fprintf(fid, '%u\n', b.neighbour-1);
fprintf(fid, '%s\n', ')');
fprintf(fid, '\n%s\n', footer);
fclose(fid);

delete([polyMeshDir filesep 'neighbour*'])
if compress
    gzip('neighbour',polyMeshDir)
    delete('neighbour')
else
    movefile('neighbour',polyMeshDir)
end

%Write boundary file
fid = fopen('boundary','w');
writeHeader(fid,'polyBoundaryMesh','boundary');

fprintf(fid, '%s\n', int2str(length(b.boundary.nFaces)));
fprintf(fid, '%s\n', '(');
for n = 1:length(b.boundary.nFaces)
    fprintf(fid, '%s\n', ['    ' b.boundary.patchName{n}]);
    fprintf(fid, '%s\n', '    {');
    fprintf(fid, '%s\n', ['        type            ' b.boundary.patchType{n} ';']);
    fprintf(fid, '%s\n', ['        nFaces          ' int2str(b.boundary.nFaces(n)) ';']);
    fprintf(fid, '%s\n', ['        startFace       ' int2str(b.boundary.startFace(n)-1) ';']);    fprintf(fid, '%s\n', '    }');
end
fprintf(fid, '%s\n\n',')');
fprintf(fid, '\n%s\n', footer);
fclose(fid);

%Boundary file is not compressed
movefile('boundary',polyMeshDir) 

function writeHeader(fid,classType,objectType,b)

str = {'/*--------------------------------*- C++ -*----------------------------------*\'
'| =========                 |                                                 |'
'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'
'|  \\    /   O peration     | Version:  2.0.1                                 |'
'|   \\  /    A nd           | Web:      www.OpenFOAM.COM                      |'
'|    \\/     M anipulation  |                                                 |'
'\*---------------------------------------------------------------------------*/'
' '
'FoamFile'
'{'
'    version         2.0;'
'    format          ascii;'
'    class           classType;'
'    location       "constant/polyMesh";'
'    object          objectType;'
'}'
' '
'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
' '};

if nargin > 3
    noteStr = [
        '    note     "nPoints: ' num2str(b.nPoints) ...
                        ' nCells: ' num2str(b.nCells) ...
                        ' nFaces: ' num2str(b.nFaces) ...
                        ' nInternalFaces: ' num2str(b.nInternalFaces) ...
                        '";'];
    str = {str{1:13} noteStr str{14:end}};
end

str = strrep(str,'classType',classType);
str = strrep(str,'objectType',objectType);

for ns = 1:length(str)
    fprintf(fid, '%s\n', str{ns});
end  