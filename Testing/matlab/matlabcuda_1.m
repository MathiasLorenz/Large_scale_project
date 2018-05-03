function matlabcuda_1(datPath,figPath)
if nargin == 0
    datPath = '../data/';
end
addpath(genpath('./'));
DataName = 'cuda_1-';
DataFiles = finddata(datPath,DataName);

%% Performance plots

for i = 1:length(DataFiles)
	FileName = [datPath,DataFiles{i}];
    M = read3DMatrixFromFile(FileName);
    testMatrix(M,DataFiles{i}(1:end-4));
end
if nargin ~= 0
ExportFigures([],figPath);
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
CreateFigure([name,'_slize']);
subplot(121);
surf(my_sol(:, :, round(0.80*Nz)));
title('Approximated solution at 80% slize')
axis([0,Nx,0,Ny,-1,1])
subplot(122);
surf(true_sol(:, :, round(0.80*Nz)));
title('Real solution at 80% slize')
axis([0,Nx,0,Ny,-1,1])

%% Plot the error

err_mat = abs(my_sol - true_sol);
CreateFigure([name,'_error']);
surf(err_mat(:, :, round(0.8*Nz)))
title('Element error for 80% slice')
xlabel('x')
ylabel('y')
zlabel('Error')

%% Plot the error through the z dimension

err_mat = abs(my_sol - true_sol);
sz = size(err_mat, 3);
err_vec = zeros(sz);

for i = 1:sz
    err_vec(i) = max(max(abs(my_sol(:, :, i) - true_sol(:, :, i))));
end
CreateFigure([name,'_maximal_error']);
plot(err_vec);
xlabel('z slize')
ylabel('Max error')
end