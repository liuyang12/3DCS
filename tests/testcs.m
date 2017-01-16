%TESTCS test compressive sensing (CS) algorithms for sparse reconstruction
%   See also CS.
%% step 1. generate and save simulation dataset
clear; clc;
% close all;
% [1.0] parameters for simulation
sparams = []; % paramters for simulation
    sparams.rows = 128; % number of rows (m)
    sparams.cols = 128; % number of columns (n)
    sparams.samprate = 0.2; % sampling rate (gamma)
    sparams.noisesnr = 50; % signal-to-noise ratio (sigma)
    sparams.sampname = 'gi'; % binary sample 
    % sparams.sampname = 'angrybird'; % binary sample 
    % sparams.sampname = 'phantom'; % simple gray-scale sample
    % sparams.sensmethod = 'binary'; % sensing matrix (multiple patterns) is binary
    sparams.sensmethod = 'mixhadamard'; % sensing matrix (multiple patterns) is binary
    sparams.SAMPLE_BINARY  = 0; % sample is binary (1-Y, 0-N)
    sparams.savesize = 5; % saving size of the image
    % sparams.binsize = 1; % binning size of the spatial light modulator
REGENDATA = true; % flag of regenerating date (1-Y,0-N)
% REGENDATA = false; % flag of regenerating date (1-Y,0-N)
if REGENDATA
% [1] generate and save dataset
    gendata(sparams);
end

%% step 2. load and apply CS method to the dataset
% [2] load dataset
clear; clc;
% close all;
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

datadir    = '../../data';       % data directory
simdatadir = '../../data/sim/static';   % simulation data directory
outdir     = '../out';           % simulation output directory
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% [2.1] load simulation preferences
load(sprintf('%s/sim_prefs.mat',datadir));
% [2.2] load dataset
load(sprintf('%s/%s%d_samprate%.2f_snr%ddb.mat',simdatadir,...
    sparams.sampname,sparams.rows,sparams.samprate,sparams.noisesnr));
samp    = csdata.samp;    % sample
sensmat = csdata.sensmat; % sensing matrixs
meas    = csdata.meas;    % measurement vector

rows    = sparams.rows;
cols    = sparams.cols;

% [3] apply CS method
csparams = []; % parameters for CS method
    csparams.rows     = rows;
    csparams.cols     = cols;
    % csparams.srbasis  = 'haar'; % sparse representation basis
    % csparams.srbasis  = 'daubechies'; % sparse representation basis
    csparams.srbasis  = 'dct'; % sparse representation basis
    
    % csparams.csmethod = 'gpsr'; % GPSR ell_1 solver
    % csparams.csmethod = 'tv'; % TV regularization method
    % csparams.csmethod = 'gap'; % GAP ell_1 solver
    % csparams.csmethod = 'gap-tv'; % GAP-TV  solver
    csparams.csmethod = 'gap2d'; % GAP2D solver
tic;
sig_out = cs(sensmat,meas,csparams); % apply cs method
t_cs = toc;
% sig_out = abs(sig_out);   % abs output
% sig_out = max(sig_out,0); % threshold output
samp_rc = reshape(sig_out,rows,cols);
samp_nm = imnorm(samp_rc);
samp_rz = imresize(samp_nm,sparams.savesize,'nearest');
recrmse = sqrt(immse(uint8(samp_nm*255),uint8(imnorm(samp)*255))); % root mean-square-error
recpsnr = psnr(uint8(samp_nm*255),uint8(imnorm(samp)*255));        % peak signal-to-noise-ratio
recssim = ssim(uint8(samp_nm*255),uint8(imnorm(samp)*255));        % structure similarity
imwrite(samp_rz,sprintf('%s/%s%d_samprate%.2f_snr%ddb_%s_rmse%.2f_psnr%.2f_ssim%.2f_t%.1f.png',outdir,...
    sparams.sampname,sparams.rows,sparams.samprate,sparams.noisesnr,csparams.csmethod,recrmse,recpsnr,recssim,t_cs));
figure; 
imshow(samp_rz);


