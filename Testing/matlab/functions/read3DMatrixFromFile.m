function [output_matrix, error] = read3DMatrixFromFile(file)
% Function to read in 3D matrix from file.
% The file must in the first line have 3 integers with Nx Ny Nz specified.
% The next line contains space seperated values of the matrix.

% Open file
fid = fopen(file, 'r');
if fid == -1
    error('Cannot open file.');
end

% Read in header
C = textscan(fid, '%d %d %d', 1);
C = double(cell2mat(C));
Nx = C(1); Ny = C(2); Nz = C(3);
total_elem = Nx*Ny*Nz;

% Read the matrix
C = textscan(fid, '%f ', total_elem);
C = cell2mat(C);

output_matrix = reshape(C, Nx, Ny, Nz);

% Close file
fclose(fid);