function matlabspeedup(DataName,BaseName,datPath,figPath)
%MATLABSPEEDUP
%  This function is designed to output a number of performance plots
%  (Memory-MFLOP plot). The inputs for the function is listed bellow.
%  __________________________________________________________________
%  DataName:
%       Base name of the files containing the data. Should be the same as
%       the jobname for LSF based jobs.
%       The data is assumed to have the name: 
%           DataName-Method.dat
%           DataName-N-Method.dat
%       Where Method is the part which will be used for the label. N defines how 
%       many cores were used for the computation. Use N if relevant,
%       otherwise ignore it.
%
%  datPath:
%       Path to where the data is located.
%
%  figPath:
%       Path to where figures should be exported to. If this is empty or
%       not defined then no figures are exported.
%
%  See also MEMORYVFLOPS.

close all
if ~exist('DataName','var') || isempty(DataName)
    DataName = 'coresu';
end
if ~exist('BaseName','var') || isempty(BaseName)
    BaseName = 'omp3d';
end


if ~exist('datPath','var') || isempty(datPath)
    datPath = '../data/';
end
addpath(genpath('./'));
DataFiles = finddata(datPath,DataName);

%% Performance plots
FileName = cell(length(DataFiles),1);
FuncName = cell(length(DataFiles),1);
for i = 1:length(DataFiles)
	FileName{i} = [datPath,DataFiles{i}];
    FuncName{i} = getFuncName(FileName{i});
end
corePerformance(FileName,BaseName,DataName,'Speedup',FuncName);
if exist('figPath','var') && ~isempty(figPath)
    ExportFigures([],figPath);
end
end

function FuncName = getFuncName(FileName)
k = strfind(FileName,'-');

if size(k,2) > 1
    Cores = FileName(k(1)+1 :k(2)-1);
    FuncName = FileName(k(2)+1:end-4);
    
    FuncName = [FuncName,' (n: ',Cores,')'];
else
	FuncName = FileName(k(1)+1:end-4);
end

spaces = strfind(FuncName,'_');
FuncName(spaces) = ' ';



end