function [Data] = loadTXT(filename)
fileID = fopen(filename);
x = textscan(fileID,'%s %f %s %f %s %f');
fclose(fileID);


Data = [x{2}, x{4}, x{6}];

