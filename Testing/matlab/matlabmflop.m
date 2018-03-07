function matlabmflop(datPath,figPath)


%% OpenMP in 2D
FuncName = {'omp2d'};
FileName = [datPath,'mflop',FuncName{1},'.dat'];
disp(FileName);
PlotName = 'omp2d_performance';
PlotTitle= 'Performance OpenMP Jacobi in 2D';
memoryVflops(FileName,PlotName,PlotTitle,FuncName);

%% OpenMP in 3D
FuncName = {'omp3d'};
FileName = [datPath,'mflop',FuncName{1},'.dat'];
PlotName = 'omp3d_performance';
PlotTitle= 'Performance OpenMP Jacobi in 3D';
memoryVflops(FileName,PlotName,PlotTitle,FuncName);

%% Compare the different methods
FileNames = {[datPath,'mflopomp2d.dat'],[datPath,'mflopomp3d.dat']};
PlotName = 'jseq_restrict';
PlotTitle= 'Restrict Vs no Restrict for Sequential Jacobi';
FuncName = {'Restrict','No Restrict'};
memoryVflopsCompare(FileNames,PlotName,PlotTitle,FuncName);

ExportFigures([],figPath);