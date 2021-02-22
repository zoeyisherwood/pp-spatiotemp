function [output] = calc_spatialslope(input,saveFlag,plotFlag,plotDir,plotName)
%
% principle code written by Colin Clifford.
% batch script put together by Zoey Isherwood.
%
% Data spit out from this function:
%
%   1) output.spatialSlope_quadFit
%   2) output.spatialSlope_logFit;
%   3) output.R2_quadFit;
%   4) output.R2_logFit;
%   5) output.chosenCutOff
%   6) output.maxPointsTested
%
% if you use this code, please cite:
% 
% Isherwood, Z. J., Clifford, C. W., Schira, M. M., Roberts, M. M. & Spehar, B. (2021) 
%     Nice and slow: Measuring sensitivity and visual preference toward naturalistic stimuli 
%     varying in their amplitude spectra in space and time. Vision Research 181, 47-60, 
%     doi:10.1016/j.visres.2021.01.001.
% 
% Log:
%
% YYYYMMDD
% 20181031 - Version 2 created.
%
% scripts needed to run this code (at end of script):
%
%   one_over_x.m by Colin Clifford
%   rot_avg.m by Bruno Olshausen
%
% before fitting the data, outliers are removed if it's Cook's Distance >
% n/4
% 
% This script will also allow non-square input... However, non-square (and
% also input that isn't a factor of 2^n) takes longer to process since
% Fourier transforms will take longer.
%
% This script will also check if your input is coloured and covert to
% grayscale if it is.
%
% To do for future versions:
%
% -Clean up plotting section of the code. A lot of redundancies.

%% input vars--------------------------------------------------------------

scriptStartTime = datestr(now, 30);

%% input checks------------------------------------------------------------

% is the input grayscale?

if ndims(input) == 3
    
    %your input is coloured... need to convert to grayscale
    
    input = rgb2gray(input);
    
end

% will the input be an integer when divided by 2? if not, remove a pixel in
% the x and y direction.

% x direction

if mod(size(input,1),2) == 1 % will be 0 if even number, 1 if odd.
    
    %if size(f,1) is odd, then when it's divided by 2 later in the
    
    %script it will no longer be an integer... to ameliorate this, we
    
    %will scrape off one pixel in both x and y directions.
    
    input = input(1:end-1,:);
    
end

% y direction

if mod(size(input,2),2) == 1 % will be 0 if even number, 1 if odd.
    
    input = input(:,1:end-1);
    
end

if size(input,1) ~= size(input,2)
    
    %     need to crop. will fix this in future versions
    
    s_original = size(input);
    s_trim = min(pow2(fix(log2(s_original(1:2)))));
    input = input(fix(s_original(1)/2) - s_trim/2+1:fix(s_original(1)/2) + s_trim/2, ...
        fix(s_original(2)/2)-s_trim/2+1:fix(s_original(2)/2)+s_trim/2,:); %taken from one_over_f.m
    
end

%% calculate spatial slope-------------------------------------------------

% compute its 2-D spectrum

spec = fftn(input);

% convert from complex to amplitude spectrum

ampspec = abs(spec);

% collapse across spatial dimensions

x = mean(ampspec,1);
xx = mean(x,2);

xsize = size(input,1);
ysize = size(input,2);
imf=fftshift(fft2(input));

impf=abs(imf);

Pf = rotavg(impf);

x_vec = 1:xsize/2;
y_vec = Pf(2:(ysize/2)+1);
y_vec = y_vec';

%% figure out what cutoff gives the best fit-------------------------------

aic_quadratic = [];
aic_loglog = [];

r2_comparison_quad = [];
r2_comparison_loglog = [];

%% choose cut off based on outliers

% Use regstats to calculate Cook's Distance
stats = regstats(y_vec,x_vec,'linear');
% if Cook's Distance > n/4 is a typical threshold that is used to suggest
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

spatialSlope_quadFit = params_0(2);

probs = one_over_x(params_0,x_vec);

% loglog

A_loglog = log(x_vec(startingValue:end));
B_loglog = log(y_vec(startingValue:end));
b_loglog = polyfit(A_loglog, B_loglog, 1);

spatialSlope_logFit = b_loglog(1);


% Calculate R2-----------------------------------------------------

% quad fit

A_quadratic = x_vec(startingValue:end);
B_quadratic = y_vec(startingValue:end);
f_quadratic = probs(startingValue:end);
Bbar_quadratic = mean(B_quadratic);
SStot_quadratic = sum((B_quadratic - Bbar_quadratic).^2);
SSres_quadratic = sum((B_quadratic - f_quadratic).^2);
R2_quadratic = 1 - SSres_quadratic/SStot_quadratic;

% linear fit (data fit on loglog axis, transformed to linear axis)

f_linear = polyval(b_loglog, A_loglog);
f_linear = exp(f_linear);
Bbar_linear = mean(B_quadratic);
SStot_linear = sum((B_quadratic - Bbar_linear).^2);
SSres_linear = sum((B_quadratic - f_linear).^2);
R2_linearSpace = 1 - SSres_linear/SStot_linear;

%% loglog fit------------------------------------------------------

A_loglog = log(x_vec);
B_loglog = log(y_vec);
b_loglog = polyfit(A_loglog, B_loglog, 1);

spatialSlope_logFit = b_loglog(1);

f_loglog = polyval(b_loglog, A_loglog);
Bbar_loglog = mean(B_loglog);
SStot_loglog = sum((B_loglog - Bbar_loglog).^2);

SSres_loglog = sum((B_loglog - f_loglog).^2);
R2_loglog = 1 - SSres_loglog/SStot_loglog;

if ~isvar('plotDir')
    
    %save in the current directory in a directory called plots
    
    dirForSaving = ['spatSlopeAnalysis'];
    
    
    if ~isdir(dirForSaving)
        
        mkdir(dirForSaving)
        
    end
    
elseif isempty(plotDir)
    
    %save in the current directory in a directory called plots
    
    dirForSaving = ['spatSlopeAnalysis'];
    
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
    
    legend(['slope ' num2str(spatialSlope_logFit)],['R^2 ' ...
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
    
    legend(['slope ' num2str(spatialSlope_quadFit)],['R^2 ' num2str(R2_quadratic)]);
    
    % linear fit on loglog axis--------------------------------------------
    
    A_linearFitLogAxis = log(x_vec(startingValue:end));
    B_linearFitLogAxis = log(y_vec(startingValue:end));
    f_linearFitLogAxis = log(probs(startingValue:end));
    Bbar_linearFitLogAxis = mean(B_linearFitLogAxis);
    SStot_linearFitLogAxis = sum((B_linearFitLogAxis - Bbar_linearFitLogAxis).^2);
    SSres_linearFitLogAxis = sum((B_linearFitLogAxis - f_linearFitLogAxis).^2);
    R2_linearFitLogAxis = 1 - SSres_linearFitLogAxis/SStot_linearFitLogAxis;
    
    subplot 133
    
    plot(A_linearFitLogAxis,B_linearFitLogAxis,'bo');
    
    hold on;
    
    plot(A_linearFitLogAxis,f_linearFitLogAxis,'r-'); % fit of a*(x^b)
    
    hold off;
    
    title(['linear fit on logged axis']);
    xlabel('log SF');
    ylabel('log amp');
    
    legend(['slope ' num2str(spatialSlope_quadFit)], ...
        ['R^2 ' num2str(R2_linearFitLogAxis)]);
    
    % save figures---------------------------------------------------------
    
    if ~isvar('plotName')
        
        plotName = [dirForSaving '/' 'spatialSlopePlot_' scriptStartTime ...
            '_cutOff' num2str(startingValue) '.png'];
        
    elseif  isempty(plotName)
        
        plotName = [dirForSaving '/' 'spatialSlopePlot_' scriptStartTime ...
            '_cutOff' num2str(startingValue) '.png'];
        
    end
    
end

% prepare everything for output...

output = struct;

output.spatialSlope_quadFit = spatialSlope_quadFit;
output.spatialSlope_logFit = spatialSlope_logFit;
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

% rotavg.m - function to compute rotational average of (square) array
% by Bruno Olshausen
%
% function f = rotavg(array)
%
% array can be of dimensions N x N x M, in which case f is of
% dimension NxM. N should be even.

function f = rotavg(array)

[N,N,M]=size(array);

[X,Y]=meshgrid(-N/2:N/2-1,-N/2:N/2-1);

[theta,rho]=cart2pol(X,Y);

rho=round(rho);
i=cell(N/2+1,1);
for r=0:N/2
    i{r+1}=find(rho==r);
end

f=zeros(N/2+1,M);

for m=1:M
    
    a=array(:,:,m);
    for r=0:N/2
        f(r+1,m)=mean(a(i{r+1}));
    end
    
end

end
