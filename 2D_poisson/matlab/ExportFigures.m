function ExportFigures(FigureList,ExportLocation,varargin)
%EXPORTFIGURES
%  This function is an automation of exporting figures. This function
%  have been modified to better fit the current MPI/CUDA project
%  __________________________________________________________________
%  EXPORTFIGURES()
%       Exports all figures currently in the Graphical Root in
%       .eps format to the standard location ../results/
%
%  EXPORTFIGURES(FigureList)
%       Exports all figures in the specified list in
%       .eps format to the standard location: ../results/
%       The list should be a list of figure variables.
%
%  EXPORTFIGURES(FigureList,ExportLocation)
%       Exports all figures in the specified list in
%       .eps format to the specified location ExportLocation
%       The list should be a list of figure variables.
%       The location should be described in a string.
%
%  EXPORTFIGURES(...,'dpi',Integer)
%       Additionally sets the dpi for the exported image.
%
%  EXPORTFIGURES(...,'ext',String)
%       Additionally defines the extension for the output image.
%       Common extension is 'png', 'jpeg', 'pdf' and 'epsc'.
%       See <a
%       href="matlab: doc print"
%       >print documentation</a> for the full list.
%
%  EXPORTFIGURES(...,'asp',Double)
%       Sets a userdefined aspect ration for the exported figure, default
%       value is 1.5.
%
%  See also GROOT, PRINT.

%% Handeling inputs
% Handle the list of figures
if ~exist('FigureList','var') || isempty(FigureList)
    FigureList = get(groot,'Children');
end

% Handle export location
if ~exist('ExportLocation','var') || isempty(ExportLocation)
    ExportLocation = '../results/';
elseif ~strcmp(ExportLocation(end),'/')
    ExportLocation = [ExportLocation,'/'];
end
if ~exist(ExportLocation,'dir')
    fprintf('\nDirectory not found, created at: %s\n',ExportLocation)
    mkdir(ExportLocation)
end

%% Handle additional inputs through varargin
if nargin > 2
    i=1;
    disp(' ')
    disp('Exporting with userdefined settings:')
    while i < length(varargin)
        switch lower(varargin{i})
            case 'dpi'
                dpi = ['-r',int2str(varargin{i+1})];
                fprintf(' - DPI: %d\n',varargin{i+1})
            case 'ext'
                ext = varargin{i+1};
                fprintf(' - EXT: %s\n',ext)
            case 'asp'
                asp = varargin{i+1};
                fprintf(' - ASPECT: %5.2f\n',asp)
            otherwise
                fprintf('Unknown option: %s\n',varargin{i});
        end
        i=i+2;
    end
end

%% Ensure existence of variables
if isempty(FigureList)
    return
end
if ~exist('dpi','var')
    dpi = '-r100';
end
if ~exist('ext','var')
    ext = 'epsc';
end
if ~exist('asp','var')
    asp = 1.5;
end


%% Exporting the figures
fprintf('Exported Figures: (Location: %s)\n',ExportLocation)
for i=1:length(FigureList)
    % Handle Name
    if isempty(FigureList(i).Name)
        FigureList(i).Name = sprintf('Figure%d',FigureList(i).Number);
    end
    
    % Handle docking, position and size
    Dock = false;
    if strcmp(FigureList(i).WindowStyle,'docked')
        Dock = true;
        FigureList(i).WindowStyle = 'normal';
    end
    FigureList(i).Position = [300,100,asp*600,600];
	drawnow;
    
    
    fprintf(' - %s.%s\n',FigureList(i).Name,ext);
    print(FigureList(i),['-d',ext],...
        [ExportLocation,FigureList(i).Name],dpi);
    
    if Dock
        FigureList(i).WindowStyle = 'docked';
    end
end
end
