function memoryVflops(FileName,PlotName,PlotTitle,FuncName,CPU)

%% Handle input
if ~exist('CPU','var')
    CPU = 'XeonE5-2660';
end

%% Get and shape the data.
[data] = loadTXT(FileName);
memory = reshape(data(:,1),length(FuncName),...
    size(data,1)/length(FuncName))';
flops  = reshape(data(:,2),length(FuncName),...
    size(data,1)/length(FuncName))';


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
F.Position = [300,100,1200,600];
title(PlotTitle,'FontSize',14)

%% Do the plots
hold on
for i = 1:size(memory,2)
    plot(log2(memory(:,i)),flops(:,i),'-d')
end
hold off
legend(FuncName,'FontSize',12)

%% Handle the axes and the ticks
n = max(2 , floor(min(log2(memory(:,1)))));
xlabel('Memory footprint')
ylabel('Mflop/s')
grid on
A = gca;
A.XLim(1) = n;
A.YLim(1) = 0;
A.XTick = floor(n);
while A.XTick(end) < A.XLim(2)
    n = n+1;
    A.XTick = [A.XTick, n];
end
for n = 1:length(A.XTick)
    if floor(2^A.XTick(n)) < 1024
        A.XTickLabel{n} = sprintf('%d kB',floor(2^A.XTick(n)));
    else
        A.XTickLabel{n} = sprintf('%d MB',floor(2^A.XTick(n)/1024));
    end
end

%% Insert Cache lines
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
