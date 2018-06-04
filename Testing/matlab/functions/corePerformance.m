function corePerformance(FileName,BaseName,PlotName,PlotTitle,FuncName)

%% Handle input
if ~exist('FuncName','var')
    FuncName = 'Dummy';
end

%% Get and shape the data.
[~,I] = max(strcmp(FuncName,BaseName));

flops = zeros(length(FileName),1);
names = cell (length(FileName),1);

[data] = loadTXT(FileName{I});
flops(1) = data(:,2);
names{1} = BaseName;
k = 2;
for i = 1:length(FileName)
    if I ~= i
        [data] = loadTXT(FileName{i});
        flops(k) = data(:,2)';
        names{k} = FuncName{i};
        k = k+1;
    end
end

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
P = bar(flops./flops(1));
hold off

%% Handle the axes and the ticks
%n = max(2 , floor(min(log2(memory(:,1)))));
xlabel('Memory footprint')
ylabel('Mflop/s')
grid on
A = gca;
A.YLim(1) = 0;
A.YLim(2) = max(max(flops./flops(1)),A.YLim(2))*1.1;
A.XTick = 1:length(names);
A.XTickLabel = names;

%% Insert Cache lines
%{
switch CPU
    case 'XeonE5-2660'
        L1 = 32;
        L2 = 256;
        L3 = 25600;
end
L1p = log2(L1);
L2p = log2(L1+L2);
L3p = log2(L1+L2+L3);
Yp  = A.YLim(2)-0.011*norm(A.YLim);
Style = {'LineStyle','--','Color','black'};
hold on
% Lines
line([L1p,L1p],A.YLim,Style{:})
line([L2p,L2p],A.YLim,Style{:})
line([L3p,L3p],A.YLim,Style{:})

% Texts
text(A.XLim(1)+0.4*(L1p-A.XLim(1)),Yp,'L1 Cache')
text(L1p+0.4*(L2p-L1p),Yp,'L2 Cache')
text(L2p+0.4*(L3p-L2p),Yp,'L3 Cache')
text(L3p+0.4*(A.XLim(2)-L3p),Yp,'RAM')
hold off
%}