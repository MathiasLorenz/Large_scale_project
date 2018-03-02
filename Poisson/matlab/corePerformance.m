function corePerformance(FileNames,PlotName,PlotTitle)
if nargin ==0
    matlabschedule
end

%% Get and shape the data.
for i = 1:length(FileNames)
    Data = loadTXT(['output',FileNames{i},'.txt']);
    Threads(:,i)=[1:16]';
    Speedup(:,i)=Data(1,3)./Data(:,3);
end
%% Fetch the plot if it exists
F = findobj('Name',PlotName);
if isempty(F)
    F = figure();
else
    figure(F);
    clf;
end

%% Name the plot and set title and size
F.Name = PlotName;
F.Position = [300,100,600,600];
title(PlotTitle,'FontSize',14)

%% Do the plots
hold on
L = line([1,max(Threads(:))],[1,max(Threads(:))]);
plot(Threads,Speedup,'--*','LineWidth',2);
legend({'theoretic',FileNames{:}},'FontSize',12)
hold off
A = gca;
A.XLim(1) = 1;
A.YLim(1) = 1;
grid on;
L.LineWidth = 1;
L.Color = 'green';