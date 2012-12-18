function meshDir = makeCaseDir(caseDir)

%Generating OpenFOAM case directory structure for paraview to be able to
%read the generated mesh (open the generated .foam file in paraview)
%
%Johan Roenby, DHI Water & Envirnment

if exist(caseDir)~=7 
    mkdir(caseDir);
end
if exist([caseDir filesep 'constant']) ~= 7
    mkdir(caseDir,'constant');
end
meshDir = [caseDir filesep 'constant' filesep 'polyMesh'];
if exist(meshDir) ~= 7
    mkdir(caseDir,['constant' filesep 'polyMesh']);
end
if exist([caseDir filesep '0']) ~= 7
    mkdir(caseDir,'0')
end
if exist([caseDir filesep 'system']) ~= 7
    copyfile([caseDir filesep 'private' filesep 'case' filesep 'system'],[caseDir filesep 'system'])
end
ind = strfind(caseDir,filesep);
caseName = caseDir(ind(end)+1:end);
fid = fopen([caseDir filesep caseName '.foam'],'w');
fclose(fid);