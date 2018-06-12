close all;
f = @(n) (n.*16.*1e9./(3*8)).^(1/3);
plot(1:12,f(1:12))
axis square;
xlabel('Number of GPU');
ylabel('Max array size');
title('Theoretical max size of array');
grid on
L = get(gca,'Children');
L.LineWidth = 2;
set(gcf,'Name','maxarraysize')
ExportFigures([],'../figures');