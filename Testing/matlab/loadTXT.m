function [Data] = loadTXT(filename)
if exist(filename,'file')
	fileID = fopen(filename);
else
	error('LOADTXT:FILENOTFOUND','Cannot locate the file:\n  %s\n',filename);
end
x = textscan(fileID,'%s %f %s %f %s %f');
fclose(fileID);


Data = [x{2}, x{4}, x{6}];

