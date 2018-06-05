function matlabparperf(DataName,datPath,figPath)
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
    DataName = 'parperf';
end


if ~exist('datPath','var') || isempty(datPath)
    datPath = '../data/';
end
addpath(genpath('./'));
DataFiles = finddata(datPath,DataName);

%% Performance plots
FileName = cell(length(DataFiles),1);

for i = 1:length(DataFiles)
    FileName{i} = [datPath,DataFiles{i}];
end

[Names,Cores] = getFuncName(FileName);

parPerformance(FileName,DataName,'Speedup',Names,Cores);
if exist('figPath','var') && ~isempty(figPath)
    ExportFigures([],figPath);
end
end

function [Names,Cores] = getFuncName(FileNames)
Cores = zeros(length(FileNames),1);
Names = cell(length( FileNames),1);
for i = 1:length(FileNames)
    FileName = FileNames{i};
    k = strfind(FileName,'-');
    
    Name = FileName(k(2)+1:end-4);
    spaces = strfind(Name,'_');
    
    Name(spaces) = ' ';
    Core = FileName(k(1)+1 :k(2)-1);
    
    Names{i} = Name;
    Cores(i) = str2double(Core);
end
end