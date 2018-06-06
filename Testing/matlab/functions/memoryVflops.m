function memoryVflops(FileName,PlotName,PlotTitle,FuncName)

%% Handle input
if ~exist('FuncName','var')
    FuncName = 'Dummy';
end

%% Get and shape the data.
[data] = loadTXT(FileName);
memory = data(:,1)';
flops  = data(:,2)'./1000;


%% Fetch the plot if it exists
F = findobj('Name',PlotName);
if isempty(F)
    F = figure();
else
    figure(F);
end

%% Name the plot and set title and size
F.Name = PlotName;
title(PlotTitle,'FontSize',14)

%% Do the plots
hold on
P = plot(log2(memory),flops,'-d');
hold off

%% Handle the legend
CurrentLegend = get(findobj(F,'Type','Legend'),'String');
if isempty(CurrentLegend)
    CurrentLegend = FuncName;
else
    CurrentLegend{end} = FuncName;
end
legend(CurrentLegend);

%% Handle the axes and the ticks
n = max(2 , floor(min(log2(memory(:,1)))));
xlabel('Memory footprint')
ylabel('Gflop/s')
grid on
A = gca;
A.XLim(1) = n;
A.YLim(1) = 0;
A.YLim(2) = max(max(flops)*1.2,A.YLim(2))*1.1;
A.YScale  = 'log';
A.XTick = floor(n);
while A.XTick(end) < A.XLim(2)
    n = n+1;
    A.XTick = [A.XTick, n];
end
for n = 1:length(A.XTick)
    if floor(2^A.XTick(n)) < 1024
        A.XTickLabel{n} = sprintf('%d kB',floor(2^A.XTick(n)));
    elseif (floor(2^A.XTick(n)) >= 1024) && (floor(2^A.XTick(n)) < 1024^2)
        A.XTickLabel{n} = sprintf('%d MB',floor(2^A.XTick(n)/1024));
    else
        A.XTickLabel{n} = sprintf('%d GB',floor(2^A.XTick(n)/1024^2));
    end
end
