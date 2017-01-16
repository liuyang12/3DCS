function [  ] = gendata3d( sparams )
%GENDATA3D generate and save three dimensioanal (3D) simulation dataset
%
% [0.0] add *utils*, *algorithms*, and *packages* to path 
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

% [0] directories configuration
datadir    = '../../data';                   % data directory
simdatadir = '../../data/sim/dynamic';       % simulation data directory
sampdir    = '../scenes/dynamic/noise_free'; % static sample directory
if ~exist(simdatadir,'dir')
    mkdir(simdatadir);
end

% [1] load or generate data for simulation
nyqnum  = sparams.rows*sparams.cols; % number of Nyquist samping (N)
sampnum = round(sparams.samprate*nyqnum); % sampling number (M)
nframe  = sparams.nframe; % number of frames for simulation (F)
width   = sparams.width; % video width (W)
height  = sparams.height; % video height (H)
format  = sparams.format; % video format ('420' for YUV 4:2:0 default)

% [1.1] load video for simulation
samppath = sprintf('%s/%s.yuv',sampdir,sparams.sampname);
mov = yuv2mov(samppath,width,height,format); % read yuv format video
nmov = length(mov);
framespath = sprintf('%s/%s%d',sampdir,sparams.sampname,sparams.rows);
if ~exist(framespath,'dir')
    mkdir(framespath);
end
if nmov<nframe % number of frames in the video less than the requested 
               % number of frames
    error('not sufficient number of frames in %s.\n',samppath);
end
for imov = 1:nframe
    cursamp = mov(imov).cdata;
    if size(cursamp,3)>1
        cursamp = rgb2gray(cursamp);
    end
    if sparams.SAMPLE_BINARY % sample is binary
        level  = graythresh(cursamp);
        cursamp = im2bw(cursamp,level);
    end
    cursamp = imnorm(double(cursamp)); % normalization to fit grayscale sample
    if width>height % crop the middle square
        left  = floor((width-height)/2)+1;
        right = floor((width+height)/2);
        cursamp_sq = cursamp(:,left:right);
    else
        left  = floor((height-width)/2)+1;
        right = floor((height+width)/2);
        cursamp_sq = cursamp(left:right,:);
    end
    cursamp_sq = imresize(cursamp_sq,[sparams.rows sparams.cols]);
    samp(:,imov) = cursamp_sq(:); % all the samples in the video (m*n)*F
    imwrite(imresize(cursamp_sq,sparams.savesize,'nearest'), sprintf('%s/frame%04d.png', framespath, imov));
end
    
% [1.2] generate sensing matrix and get measurement vector as dataset
switch sparams.sensmethod
    case 'binary' % binary pseudo-random sensing matrix (A)
        sensmat = rand(sampnum,nyqnum,nframe)>0.5; % binary sensing matrix (A)
    case 'gaussian' % Gaussian sensing matrix (A)
        sensmat = randn(sampnum,nyqnum,nframe); 
        % sensmat = 1/sampnum*randn(sampnum,nyqnum); % Gaussian sensing matrix (A)
    case 'mixhadamard' % mixed Hadamard sensing method
        Hr = hadamard(sparams.rows); % Hadamard matrix in rows
        Hc = hadamard(sparams.cols); % Hadamard matrix in columns
        Had = kron(Hc',Hr');
        for iframe = 1:nframe
            Rad = double(rand(nyqnum,1)>0.5)*2-1; % Rademacher distribution
            allrowset = randperm(nyqnum);
            rowset = allrowset(randperm(sampnum));
            sensmat(:,:,iframe) = Had(rowset,:)*diag(Rad);
            % csdata.Rad    = Rad;
            % csdata.rowset = rowset;
        end
    otherwise
        error('unsupported sensing method.\n');
end
meas_nonoise = fm_3dcs(sensmat,samp); % forward model of 3D-CS
noisenorm  = randn(sampnum,nframe);
noise = norm(meas_nonoise,2)/norm(noisenorm,2)*10^(-sparams.noisesnr/20)*noisenorm;
meas = meas_nonoise+noise; % measurements (y)

% [1.3] save generated dataset for simulation
csdata.samp    = samp;
csdata.sensmat = sensmat;
csdata.meas    = meas;
save(sprintf('%s/%s%dby%d_samprate%.2f_snr%ddb.mat',simdatadir,...
    sparams.sampname,sparams.rows,sparams.nframe,sparams.samprate,...
    sparams.noisesnr),'csdata');
% [1.4] save sparams as simulation preferences
save(sprintf('%s/sim_prefs.mat',datadir),'sparams');

end

