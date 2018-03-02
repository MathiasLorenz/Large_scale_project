function output_matrix = read3DMatrixFromFile(file)

% Open file
fid = fopen(file, 'r');

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