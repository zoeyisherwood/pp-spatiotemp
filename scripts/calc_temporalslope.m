function [output] = calc_temporalslope(input,saveFlag,plotFlag,plotDir,plotName)
%
% principle code written by Colin Clifford.
% batch script put together by Zoey Isherwood.
%
% Data spit out from this function:
%
%   1) output.temporalSlope_quadFit
%   2) output.temporalSlope_logFit;
%   3) output.R2_quadFit;
%   4) output.R2_logFit;
%   5) output.chosenCutOff
%   6) output.maxPointsTested
%
% if you use this code, please cite:
% 
% Isherwood, Z. J., Clifford, C. W. G., Schira, M. M., Roberts, M. M. & Spehar, B. (2021) 
%     Nice and slow: Measuring sensitivity and visual preference toward naturalistic stimuli 
%     varying in their amplitude spectra in space and time. Vision Research 181, 47-60, 
%     doi:10.1016/j.visres.2021.01.001.
% 
% Log:
%
% YYYYMMDD
% 20200416 - Version 1 created.
%
% scripts needed to run this code (at end of script):
%
%   one_over_x.m by Colin Clifford
%   (optional: load_frames.m by Zoey Isherwood - will load folder of images
%   into a variable that can be used as 'input' in the current script).
%
% before fitting the data, outliers are removed if it's Cook's Distance >
% n/4
%
% To do for future versions:
%
% -Clean up plotting section of the code. A lot of redundancies.

%% input vars--------------------------------------------------------------

scriptStartTime = datestr(now, 30);

%% input checks------------------------------------------------------------

% too add.

xsize = size(input,1);
ysize = size(input,2);
zsize = size(input,3);

%% calculate temporal slope-------------------------------------------------

spec = fftn(input);

% convert from complex to amplitude spectrum
ampspec = abs(spec);

% collapse across x & y dimensions

x = mean(ampspec,1);
xx = mean(x,2);
xxx = squeeze(xx);

x_vec = 1:(zsize/2);
y_vec = xxx(2:((zsize/2)+1));
y_vec = y_vec';

%% figure out what cutoff gives the best fit-------------------------------

aic_quadratic = [];
aic_loglog = [];

r2_comparison_quad = [];
r2_comparison_loglog = [];

%% choose cut off based on outliers

% Use regstats to calculate Cook's Distance
stats = regstats(y_vec,x_vec,'linear');
% if Cook's Distance > n/4 is a typical treshold that is used to suggest
% the presence of an outlier
potential_outlier = stats.cookd > 4/length(x_vec);

noutliers = numel(potential_outlier(potential_outlier == 1));

x_vec = x_vec(imcomplement(potential_outlier));
y_vec = y_vec(imcomplement(potential_outlier));

startingValue = 1; %already got rid of the outliers

% recalc the slope and other variables for plotting for the
% determined cutOff------------------------------------------------

% linear

% fit hyperbolic to raw data ...
maxParam = max(y_vec(startingValue:end));
init_params = [maxParam -1]; % first guess to start iteration ...


params_0 = nlinfit(x_vec(startingValue:end),y_vec(startingValue:end), ...
    'one_over_x',init_params);

temporalSlope_quadFit = params_0(2);

probs = one_over_x(params_0,x_vec);

% loglog

A_loglog = log(x_vec(startingValue:end));
B_loglog = log(y_vec(startingValue:end));
b_loglog = polyfit(A_loglog, B_loglog, 1);

temporalSlope_logFit = b_loglog(1);


% Calculate R2-----------------------------------------------------

% quad fit

A_quadratic = x_vec(startingValue:end);
B_quadratic = y_vec(startingValue:end);
f_quadratic = probs(startingValue:end);
Bbar_quadratic = mean(B_quadratic);
TStot_quadratic = sum((B_quadratic - Bbar_quadratic).^2);
TSres_quadratic = sum((B_quadratic - f_quadratic).^2);
R2_quadratic = 1 - TSres_quadratic/TStot_quadratic;

% linear fit (data fit on loglog axis, transformed to linear axis)

f_linear = polyval(b_loglog, A_loglog);
f_linear = exp(f_linear);
Bbar_linear = mean(B_quadratic);
TStot_linear = sum((B_quadratic - Bbar_linear).^2);
TSres_linear = sum((B_quadratic - f_linear).^2);
R2_linearSpace = 1 - TSres_linear/TStot_linear;

%% loglog fit------------------------------------------------------

A_loglog = log(x_vec);
B_loglog = log(y_vec);
b_loglog = polyfit(A_loglog, B_loglog, 1);

temporalSlope_logFit = b_loglog(1);

f_loglog = polyval(b_loglog, A_loglog);
Bbar_loglog = mean(B_loglog);
TStot_loglog = sum((B_loglog - Bbar_loglog).^2);

TSres_loglog = sum((B_loglog - f_loglog).^2);
R2_loglog = 1 - TSres_loglog/TStot_loglog;

if ~isvar('plotDir')
    
    %save in the current directory in a directory called plots
    
    dirForSaving = ['temporalSlopeAnalysis'];
    
    
    if ~isdir(dirForSaving)
        
        mkdir(dirForSaving)
        
    end
    
elseif isempty(plotDir)
    
    %save in the current directory in a directory called plots
    
    dirForSaving = ['temporalSlopeAnalysis'];
    
    mkdir(dirForSaving)
    
    if ~isdir(dirForSaving)
        
        mkdir(dirForSaving)
        
    end
    
else
    
    dirForSaving = plotDir;
    
    if ~isdir(dirForSaving)
        
        mkdir(dirForSaving)
        
    end
    
end


if plotFlag == 1
    
    % loglog fit-----------------------------------------------------------
    
    figure;
    
    subplot 131
    
    plot(A_loglog,B_loglog,'bo');
    
    hold on
    
    plot([A_loglog(1), A_loglog(end)], b_loglog(2) ...
        + b_loglog(1).*[A_loglog(1),A_loglog(end)],'r-');
    
    hold off
    
    title(['log fit on log axis']);
    xlabel('log SF');
    ylabel('log amp');
    
    legend(['slope ' num2str(temporalSlope_logFit)],['R^2 ' ...
        num2str(R2_loglog)]);
    
    % linear fit on linear axis--------------------------------------------
    
    subplot 132
    
    hold on;
    
    plot(A_quadratic,B_quadratic,'bo'); % data
    hold on;
    
    plot(A_quadratic,f_quadratic,'r-'); % fit of a*(x^b)
    hold off;
    
    title(['linear fit on linear axis']);
    xlabel('SF');
    ylabel( 'amp');
    
    legend(['slope ' num2str(temporalSlope_quadFit)],['R^2 ' num2str(R2_quadratic)]);
    
    % linear fit on loglog axis--------------------------------------------
    
    A_linearFitLogAxis = log(x_vec(startingValue:end));
    B_linearFitLogAxis = log(y_vec(startingValue:end));
    f_linearFitLogAxis = log(probs(startingValue:end));
    Bbar_linearFitLogAxis = mean(B_linearFitLogAxis);
    TStot_linearFitLogAxis = sum((B_linearFitLogAxis - Bbar_linearFitLogAxis).^2);
    TSres_linearFitLogAxis = sum((B_linearFitLogAxis - f_linearFitLogAxis).^2);
    R2_linearFitLogAxis = 1 - TSres_linearFitLogAxis/TStot_linearFitLogAxis;
    
    subplot 133
    
    plot(A_linearFitLogAxis,B_linearFitLogAxis,'bo');
    
    hold on;
    
    plot(A_linearFitLogAxis,f_linearFitLogAxis,'r-'); % fit of a*(x^b)
    
    hold off;
    
    title(['linear fit on logged axis']);
    xlabel('log SF');
    ylabel('log amp');
    
    legend(['slope ' num2str(temporalSlope_quadFit)], ...
        ['R^2 ' num2str(R2_linearFitLogAxis)]);
    
    % save figures---------------------------------------------------------
    
    if ~isvar('plotName')
        
        plotName = [dirForSaving '/' 'temporalSlopePlot_' scriptStartTime ...
            '_cutOff' num2str(startingValue) '.png'];
        
    elseif  isempty(plotName)
        
        plotName = [dirForSaving '/' 'temporalSlopePlot_' scriptStartTime ...
            '_cutOff' num2str(startingValue) '.png'];
        
    end
    
end

% prepare everything for output...

output = struct;

output.spatialSlope_quadFit = temporalSlope_quadFit;
output.spatialSlope_logFit = temporalSlope_logFit;
output.R2_quadFit = R2_quadratic;
output.R2_logFit = R2_loglog;
output.chosenCutOff = startingValue;
output.noutliers = noutliers;

% save output

if saveFlag == 1
    
    dataSaveName = [dirForSaving '/' 'spatialSlopeOutput_' scriptStartTime ...
        '_cutOff' num2str(startingValue) '.mat'];
    
    save(dataSaveName,'output')
    
    saveas(gcf,plotName)
    
    close(gcf)
    
end

end

function probs = one_over_x(params,x)

% Design matrix for f(x) = a*x^b.
%
%
% % Colin Clifford 6.3.01

offSetParam = 0; %100; %fit rounding later down the track...

probs = params(1).*(x.^(params(2)))+offSetParam;

end

%--------------------------------------------------------------------------

function [frames] = load_frames(frameDir,imExt,grayscale)
%
% This function will load all the frames contained in the input directory.
% imExt refers to the extension of the images you wish to load (i.e.
% 'bmp','png')
%
% set the third function input term to 1 if you want to load your frames in
% grayscale
%
% Dependent scripts:
%     find_file_path - this script will get all the file names/full path
%     for loading
%
% Log
%   YYYMMDD
%   20181024 - Initialised. zji.


% first find all the files to load

[fileLocations] = find_file_path(frameDir,imExt);

% get x y z size by loading in one frame first

workingVar = imread(fileLocations{1}.fullFilePath);


% now load all the files into a matrix

if grayscale == 0 && ndims(workingVar == 3)
    
    movieSize = [size(workingVar,1) size(workingVar,2), ...
        size(workingVar,3), numel(fileLocations)];
    
else
    
    movieSize = [size(workingVar,1) size(workingVar,2) numel(fileLocations)];
    
end

if isa(workingVar,'uint8') == 1
    
    frames = uint8(zeros(movieSize));
    
else
    
    frames = zeros(movieSize);
    
end

for numFrame = 1:numel(fileLocations)
    
    if grayscale == 1 && ndims(workingVar) == 3
        
        frames(:,:,numFrame) = rgb2gray ...
            (imread(fileLocations{numFrame}.fullFilePath));
        
    else
        
        frames(:,:,numFrame) = imread(fileLocations{numFrame}.fullFilePath);
        
    end
    
end

end



