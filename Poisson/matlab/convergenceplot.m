function convergenceplot(FileBase,N,PlotTitle,PlotName)
addpath('../data/')
Error = zeros(length(N),1);
for i = 1:length(N)
    File = [FileBase,num2str(N(i)), '.txt'];
    mat = importdata(File);
    if isstruct(mat)
        dataMat = mat.data;
    else
        dataMat = mat;
    end
    [s1, ~] = size(dataMat);
    x = linspace(-1, 1, s1);
    y = x;
    [X, Y] = meshgrid(x,y);

    realsol = sin(pi.*X).*sin(pi.*Y);

    Error(i) = max(abs(dataMat(:)-realsol(:)));
end

F = findobj('Name',PlotName);
if isempty(F)
    F = figure();
else
    figure(F);
    clf;
end

%% Name the plot and set title and size
F.Name = PlotName;
F.Position = [300,100,1200,600];
title(PlotTitle,'FontSize',14)

%% Do the plots
hold on
plot(N,Error,'--*');
plot(N,1./N.^2);
hold off
legend({'Error','O(N^2)'},'FontSize',12);

%% Style plot
A = gca;
xlabel('Number of grid points')
ylabel('Error')
grid 'on';
A.XScale = 'log';
A.YScale = 'log';





