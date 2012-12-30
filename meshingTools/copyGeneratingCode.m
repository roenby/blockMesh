function copyGeneratingCode(meshDir,toolboxDir,mfile)

%Copy m-files used to generate a polyMesh into the case folder 
%case/constant/polyMesh/history/<timestamp>
%
%Johan Roenby, DHI Water & Envirnment

histDir = [meshDir filesep 'history'];
if exist(histDir) ~= 7
    mkdir(histDir);
end
time = datestr(now,30);
copyDir = [histDir filesep time];
mkdir(copyDir)
copyfile(toolboxDir,[copyDir filesep 'meshingTools'])
ind = strfind(mfile,filesep);
mfilnam = mfile(ind(end)+1:end);
copyfile([mfile '.m'],[copyDir filesep mfilnam '.m'] )