% Script to load in 3D solution, calculate the error and plotcl
clc, clear, close all
addpath 'functions'

fname = '../../Poisson/200grid.txt';

my_sol = read3DMatrixFromFile(fname);
[Nx,Ny,Nz] = size(my_sol);
gvx = linspace(-1, 1, Nx);
gvy = linspace(-1, 1, Ny);
gvz = linspace(-1, 1, Nz);
[X, Y, Z] = meshgrid(gvx, gvy, gvz);

true_sol = sin(pi*X).*sin(pi*Y).*sin(pi*Z);

%% Calculate max error
[err, err_idx] = max(abs(my_sol(:) - true_sol(:)));
disp(err)


%% Plot the solution with true solution
CreateFigure('25%');
subplot(121);
surf(my_sol(:, :, round(0.25*Nz) ));
axis([0,Nx,0,Ny,-1,1])
subplot(122);
surf(true_sol(:, :, round(0.25*Nz)));
axis([0,Nx,0,Ny,-1,1])

CreateFigure('50%');
subplot(121);
surf(my_sol(:, :, round(0.50*Nz)));
axis([0,Nx,0,Ny,-1,1])
subplot(122);
surf(true_sol(:, :, round(0.50*Nz)));
axis([0,Nx,0,Ny,-1,1])

CreateFigure('80%');
subplot(121);
surf(my_sol(:, :, round(0.80*Nz)));
axis([0,Nx,0,Ny,-1,1])
subplot(122);
surf(true_sol(:, :, round(0.80*Nz)));
axis([0,Nx,0,Ny,-1,1])

CreateFigure('100%');
subplot(121);
surf(my_sol(:, :, round(Nz)));
axis([0,Nx,0,Ny,-1,1])
subplot(122);
surf(true_sol(:, :, round(Nz)));
axis([0,Nx,0,Ny,-1,1])

%% Plot the error

err_mat = abs(my_sol - true_sol);
CreateFigure('Error at 70%');
surf(err_mat(:, :, round(0.34*Nz)))

%% 
close all


%% Plot the error through the z dimension

err_mat = abs(my_sol - true_sol);
sz = size(err_mat, 3);
err_vec = zeros(sz);

for i = 1:sz
    err_vec(i) = max(max(abs(my_sol(:, :, i) - true_sol(:, :, i))));
end
CreateFigure('Error at through z');
plot(err_vec);