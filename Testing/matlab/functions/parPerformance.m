function parPerformance(FileName,PlotName,PlotTitle,Names,Cores)

%% Get and shape the data.

flops = zeros(length(FileName),1);

for i = 1:length(FileName)
    [data] = loadTXT(FileName{i});
    flops(i) = data(:,2)';
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
N = unique(Names);
C = unique(Cores);
Data = zeros(length(C),length(N));

for i = 1:length(Names)
    for n = 1:length(N)
        if strcmp(Names{i},N{n})
            Data(Cores(i)-1,n) = flops(i);
            
            fprintf('Name: %s, C: %d, flops: %f\n',Names{i},Cores(i),flops(i))
        end
    end
end
hold on
for n = 1:length(N)
    plot(C,Data(:,n)./Data(1,n))
end
plot(2:C(end)+2,1:C(end)+1,'LineStyle','--')
hold off

%% Handle the axes and the ticks
%n = max(2 , floor(min(log2(memory(:,1)))));
xlabel('Number of cores')
ylabel('Relative speedup')
grid on
A = gca;
A.YLim = [0,C(end)];
axis square

legend(unique(Names),'Optimal')
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