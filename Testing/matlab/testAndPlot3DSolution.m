% Script to load in 3D solution, calculate the error and plotcl
clc, clear, close all

fname = '../../Poisson/tmp2.dat';

my_sol = read3DMatrixFromFile(fname);

gvx = linspace(-1, 1, size(my_sol,1));
gvy = linspace(-1, 1, size(my_sol,2));
gvz = linspace(-1, 1, size(my_sol,3));
[X, Y, Z] = meshgrid(gvx, gvy, gvz);

true_sol = sin(pi*X).*sin(pi*Y).*sin(pi*Z);

%% Calculate max error
[err, err_idx] = max(abs(my_sol(:) - true_sol(:)));
disp(err)


%% Plot the solution with true solution
close all
figure;
surf(my_sol(:, :, 25));
figure;
surf(true_sol(:, :, 25));

figure;
surf(my_sol(:, :, 75));
figure;
surf(true_sol(:, :, 75));

figure;
surf(my_sol(:, :, 155));
figure;
surf(true_sol(:, :, 155));

%% Plot the error
close all

err_mat = abs(my_sol - true_sol);

surf(err_mat(:, :, 40))

%% Plot the error through the z dimension
close all

err_mat = abs(my_sol - true_sol);
sz = size(err_mat, 3);
err_vec = zeros(sz);

for i = 1:sz
    err_vec(i) = max(max(abs(my_sol(:, :, i) - true_sol(:, :, i))));
end

plot(err_vec);