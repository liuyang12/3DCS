function [  ] = gendata( sparams )
%GENDATA generate and save simulation dataset
%
% [0.0] add *utils*, *algorithms*, and *packages* to path 
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

% [0] directories configuration
datadir    = '../../data';       % data directory
simdatadir = '../../data/sim/static';   % simulation data directory
sampdir    = '../scenes/static'; % static sample directory
if ~exist(simdatadir,'dir')
    mkdir(simdatadir);
end

% [1] load or generate data for simulation
nyqnum  = sparams.rows*sparams.cols; % number of Nyquist samping (N)
sampnum = round(sparams.samprate*nyqnum); % sampling number (M)

% [1.1] load sample for simulation
samppath = sprintf('%s/%s',sampdir,sparams.sampname);
samp = imreadallfmt(samppath); % read sample image
if size(samp,3)>1
    samp = rgb2gray(samp);
end
if sparams.SAMPLE_BINARY % sample is binary
    level = graythresh(samp);
    samp = im2bw(samp, level);
end
samp = imnorm(double(samp)); % normalization to fit grayscale sample
samp = imresize(samp, [sparams.rows sparams.cols]);
imwrite(imresize(samp,sparams.savesize,'nearest'), sprintf('%s/%s%d.png', sampdir, sparams.sampname, sparams.rows));

% [1.2] generate sensing matrix and get measurement vector as dataset
switch sparams.sensmethod
    case 'binary' % binary pseudo-random sensing matrix (A)
        sensmat = rand(sampnum,nyqnum)>0.5; % binary sensing matrix (A)
    case 'gaussian' % Gaussian sensing matrix (A)
        sensmat = randn(sampnum,nyqnum); 
        % sensmat = 1/sampnum*randn(sampnum,nyqnum); % Gaussian sensing matrix (A)
    case 'mixhadamard' % mixed Hadamard sensing method
        Hr = hadamard(sparams.rows); % Hadamard matrix in rows
        Hc = hadamard(sparams.cols); % Hadamard matrix in columns
        Had = kron(Hc',Hr');
        Rad = double(rand(nyqnum,1)>0.5)*2-1; % Rademacher distribution
        allrowset = randperm(nyqnum);
        rowset = allrowset(randperm(sampnum));
        sensmat = Had(rowset,:)*diag(Rad);
        % csdata.Rad    = Rad;
        % csdata.rowset = rowset;
    otherwise
        error('unsupported sensing method.\n');
end
meas_nonoise = sensmat*samp(:);
noisenorm  = randn(sampnum,1);
noise = norm(meas_nonoise,2)/norm(noisenorm,2)*10^(-sparams.noisesnr/20)*noisenorm;
meas = meas_nonoise+noise; % measurements (y)

% [1.3] save generated dataset for simulation
csdata.samp    = samp;
csdata.sensmat = sensmat;
csdata.meas    = meas;
save(sprintf('%s/%s%d_samprate%.2f_snr%ddb.mat',simdatadir,sparams.sampname,sparams.rows,sparams.samprate,sparams.noisesnr),'csdata');
% [1.4] save sparams as simulation preferences
save(sprintf('%s/sim_prefs.mat',datadir),'sparams');

end

