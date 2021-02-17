function f = make_fractal_3D(xexpo,texpo,imsize,frames,contrast,aperture)
%
% Generates fractal textures with given exponents of the spatial & temporal
% amplitude spectrum & RMS contrast
%
%       xexpo - Spatial Exponent (alpha) of 1/f^alpha amplitude spectrum
%       (SS - Spatial Slope).
%
%       texpo - Temporal Exponent (alpha) of the 1/f^alpha amplitude
%       spectrum. (TS - Temporal Slope).
%
%       imsize - pixel size of your image. x an y dim (the image will be
%       square). recommend values of 2^N. i.e. 64,128,256 etc. (works
%       better with FFTs)
%
%       frames - number of frames. also recommend values of 2^N (works
%       better with FFTs).
%
%       contrast - RMS contrast. Input value between 0 & 1 (i.e 0.3 for 30%
%       RMS contrast)
%
%       aperture - 0: no aperture, 1: aperture. If set to 1, your stimulus
%       will have a spatial and temporal aperture. spatial - circular
%       gaussian. temporal - cosine ramping stimulus onset and offset.
%
%  % example:
%       make_fractal_3D(1.25,2.25,128,128,0.3,0)
%
% Code initially written by Colin Clifford. Edited by Zoey Isherwood.
%
% CC - 2.5.13 & 15.7.14 - modified from make_fractal_CC to give movie
% rather than static image
%
% ZJI - 19.7.19 - added flag to make aperture optional. set to 1 if you
% want cosine ramping and a circular aperture around the stimulus. updated
% movie function from movie2avi to VideoWriter. movie2avi only works with
% older versions of Matlab (2012 and lower).

%%
xsize = imsize;
ysize = imsize;

a = zeros(xsize,ysize,frames);
b = zeros(xsize,ysize,frames);
c = zeros(xsize,ysize,frames);
e = zeros(xsize,ysize,frames);

a = random('Normal',0,1,[xsize,xsize,frames]); % generate image of Gaussian white noise
dc = mean(mean(mean(a))); % set the mean level to zero ...
b = fftn(a-dc);     % ... and Fourier transform to get the frequency spectrum

b = fftshift(b);    % rearrange the frequency spectrum so that zero is at the centre

x0 = (xsize/2)+1; %(xsize+1)/2; %
y0 = (ysize/2)+1; %(ysize+1)/2; %
t0 = (frames/2)+1; %(frames+1)/2; %

%b(x0,y0,t0) = 0;

% h=waitbar(0,'create an amplitude spectrum with the desired "fractal" drop-off with frequency...');
disp('create an amplitude spectrum with the desired "fractal" drop-off with frequency......');

for x = 1:xsize     % create an amplitude spectrum with the desired "fractal" drop-off with frequency
    waitbar(x/xsize);
    for y = 1:ysize
        d = sqrt((x-x0)^2+(y-y0)^2);
        for t = 1:frames
            tt = sqrt((t - t0)^2);
            c(x,y,t) = (d.^(-xexpo)).*(tt.^(-texpo));
        end
    end
end
% close(h)

disp('done.')

% set DC component to zero
c(x0,y0,:) = 0;
c(:,:,t0) = 0;

b = b.*c;   % multiply the frequency spectrum of your noise by the fractal amplitude spectrum to get fractal noise

b = ifftshift(b); % rearrange the frequency spectrum so that zero is in the corner

e = ifftn(b); % inverse Fourier transform your spectrum to get a fractal noise image sequence

% normalize image sequence ...
f = real(e);
maxf = max(max(max(f)));
minf = min(min(min(f)));
ampf = max(maxf,-minf);
f = f./ampf;
ff = reshape(f,1,imsize*imsize*frames);
f = f./std(ff); % mean is 0, std is 1

cont = ones(1,frames); % initialize to peak dot contrast (defined above)

if aperture == 1
    
    % define contrast window here ...
    smooth = imsize/32;
    % define raised cosine temporal window over first & last 100ms
    temporal_smooth = 8;
    
    win_length = round(frames/temporal_smooth);     % define window length
    cont(1:win_length) = 0.5.*(1-cos(pi*(0:(win_length-1))./win_length)); % ramp up
    cont(end:-1:end-win_length+1) = cont(1:win_length); % ... and down
    
end

% h=waitbar(0,'define raised cosine temporal window over first & last 100ms...');

if aperture == 1
    disp('define raised cosine temporal window over first & last 100ms......');
else
    disp('contrast controlling stimulus...')
end

for x = 1:xsize
    waitbar(x/xsize);
    for y = 1:ysize
        
        if aperture == 1
            
            d = sqrt((x-x0)^2+(y-y0)^2);
            
            if (d > (xsize/2-1))
                spat_cont = 0;
            elseif (d > (xsize/2-1-smooth))
                spat_cont = cos(pi.*(d -(xsize/2-1)+smooth)/(2*smooth));
            else
                spat_cont = 1;
            end
            contrast_window(x,y,:) = cont.*spat_cont.*contrast;
        else
            contrast_window(x,y,:) = cont.*contrast;
            
        end
    end
end
% close(h)

disp('done.')

% implement spatio-temporal windowing here by using contrast_window r.t.
% contrast
f = min(255,max(127.5.*(1 + f.*contrast_window),0)); % clip to 0-255 with chosen contrast

g = uint8(f); % convert to unsigned 8-bit integer and write as avi ...

%% make movie file:

%--------------------------------------------------------------------------
% cc's old code:
% define colourmap for movie x = 0:(1/255):1; cm = [x;x;x]';
%
% h=waitbar(0,'Creating frames for movie file...'); disp('Creating frames
% for movie file......'); for mm = 1:frames
%     waitbar(mm/frames); M(mm).cdata = g(:,:,mm); %make movie file
%     M(mm).colormap = cm;
% end
%
% close(h)

% pathname = '..\Movies\fractmp5'; disp('executing movie2avi...')
% movie2avi(M,pathname,'compression','None','fps',60); %must remove file
% compression and set fps to 60 (the frame rate)

%--------------------------------------------------------------------------
% added by zji:

% create file name:
timeStamp = datestr(now,30);
ss_name = num2str(xexpo);
ss_name(ss_name == '.') = '-';
ts_name = num2str(xexpo);
ts_name(ts_name == '.') = '-';

fileName = [pwd '/' 'SS' ss_name '_' 'TS' ts_name '_' 'imsize' '_' num2str(imsize) ...
    '_' 'nFrames' '_' num2str(frames) '_' timeStamp '.avi' ];

% Prepare the video file...
vidObj = VideoWriter(fileName,'Uncompressed AVI');
open(vidObj);

%write frames sequentially to file....

% h=waitbar(0,'Creating frames for movie file...');
disp('Creating frames for movie file......');
for mm = 1:frames
    waitbar(mm/frames);
    currFrame = g(:,:,mm); %make movie file
    writeVideo(vidObj,currFrame);
end

close(vidObj)
% close(h)
disp('Stimulus Generated!')

% pcolor(f(:,:,33)); shading flat; colormap gray; axis equal; caxis([0
% 255]) axis off;

