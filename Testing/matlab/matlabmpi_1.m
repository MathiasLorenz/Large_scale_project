function matlabmpi_1(datPath,figPath)
if nargin == 0
    datPath = '../data/';
end
addpath(genpath('./'));
DataName = 'mpi_1-';
DataFiles = finddata(datPath,DataName);
PlotTitle= {'Performance OpenMP Jacobi in 2D',...
			'Performance OpenMP Jacobi in 3D'};

%% Performance plots

for i = 1:length(DataFiles)
	FileName = [datPath,DataFiles{i}];
    M = read3DMatrixFromFile(FileName);
    testMatrix(M,DataFiles{i}(length(DataName)+1:end));
end
end



function testMatrix(my_sol,name)

[Nx,Ny,Nz] = size(my_sol);
gvx = linspace(-1, 1, Nx);
gvy = linspace(-1, 1, Ny);
gvz = linspace(-1, 1, Nz);
[X, Y, Z] = meshgrid(gvx, gvy, gvz);

true_sol = sin(pi*X).*sin(pi*Y).*sin(pi*Z);

%% Calculate max error
err = max(abs(my_sol(:) - true_sol(:)));


%% Plot the solution with true solution
CreateFigure([name,' 25%']);
subplot(121);
surf(my_sol(:, :, round(0.25*Nz) ));
axis([0,Nx,0,Ny,-1,1])
subplot(122);
surf(true_sol(:, :, round(0.25*Nz)));
axis([0,Nx,0,Ny,-1,1])

CreateFigure([name,' 50%']);
subplot(121);
surf(my_sol(:, :, round(0.50*Nz)));
axis([0,Nx,0,Ny,-1,1])
subplot(122);
surf(true_sol(:, :, round(0.50*Nz)));
axis([0,Nx,0,Ny,-1,1])

CreateFigure([name,' 80%']);
subplot(121);
surf(my_sol(:, :, round(0.80*Nz)));
axis([0,Nx,0,Ny,-1,1])
subplot(122);
surf(true_sol(:, :, round(0.80*Nz)));
axis([0,Nx,0,Ny,-1,1])

CreateFigure([name,' 100%']);
subplot(121);
surf(my_sol(:, :, round(Nz)));
axis([0,Nx,0,Ny,-1,1])
subplot(122);
surf(true_sol(:, :, round(Nz)));
axis([0,Nx,0,Ny,-1,1])

%% Plot the error

err_mat = abs(my_sol - true_sol);
CreateFigure([name,' Error at 70%']);
surf(err_mat(:, :, round(1*Nz)))

%% Plot the error through the z dimension

err_mat = abs(my_sol - true_sol);
sz = size(err_mat, 3);
err_vec = zeros(sz);

for i = 1:sz
    err_vec(i) = max(max(abs(my_sol(:, :, i) - true_sol(:, :, i))));
end
CreateFigure([name,' Error at through z']);
plot(err_vec);
end