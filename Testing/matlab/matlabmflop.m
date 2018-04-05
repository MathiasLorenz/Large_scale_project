function matlabmflop(datPath,figPath)
addpath(genpath('./'));
FuncName = {'omp2d','omp3d'};
PlotTitle= {'Performance OpenMP Jacobi in 2D',...
			'Performance OpenMP Jacobi in 3D'};

%% Performance plots

for i = 1:length(FuncName)
	FileName = [datPath,'mflop',FuncName{i},'.dat'];
	PlotName = [FuncName{i},'_performance'];
	memoryVflops(FileName,PlotName,PlotTitle{i},FuncName{i});
end

ExportFigures([],figPath);