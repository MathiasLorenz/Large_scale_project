function matlabmflop
addpath('../results');

%% OpenMP in 2D
FuncName = {'omp2d'};
FileName = ['mflop',FuncName{1}'.txt'];
PlotName = 'omp2d_performance';
PlotTitle= 'Performance OpenMP Jacobi in 2D';
memoryVflops(FileName,PlotName,PlotTitle,FuncName);

%% OpenMP in 3D
FuncName = {'omp3d'};
FileName = 'mflopomp3d.txt';
PlotName = 'omp3d_performance';
PlotTitle= 'Performance OpenMP Jacobi in 3D';
memoryVflops(FileName,PlotName,PlotTitle,FuncName);

%% Compare the different methods
FileNames = {'mflopomp2d.txt','mflopomp3d.txt'};
PlotName = 'jseq_restrict';
PlotTitle= 'Restrict Vs no Restrict for Sequential Jacobi';
FuncName = {'Restrict','No Restrict'};
memoryVflopsCompare(FileNames,PlotName,PlotTitle,FuncName);