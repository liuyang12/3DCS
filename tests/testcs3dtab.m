%TESTCS3DTAB Test compressive sensing (CS) algorithms for three dimensional
%sparse reconstruction to construct a table of performances varing sampling
%rates.
%   See also CS3D, TEST3DCS.
%% step 1. generate and save simulation dataset
clear; clc;
% close all;
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

datadir    = '../../data';               % data directory
simdatadir = '../../data/sim/dynamic';   % simulation data directory
outdir     = '../vout';                  % simulation output directory
if ~exist(outdir,'dir')
    mkdir(outdir);
end

allsamprate = [10 20 40 60 80]*1e-2;
all3dmethod = {'tv','gap','tv3d','gap3d','gap_tv3d'};

for isamp = 1:length(allsamprate)
    % [1.0] parameters for simulation
    sparams = []; % paramters for simulation
        sparams.rows     = 64; % number of rows (m)
        sparams.cols     = 64; % number of columns (n)
        sparams.nframe   = 8;   % number of frames (F)
        % sparams.samprate = 0.4; % sampling rate (gamma)
        sparams.samprate = allsamprate(isamp); % sampling rate (gamma)
        sparams.noisesnr = 50;  % signal-to-noise ratio (sigma)
        sparams.sampname = 'foreman'; % grayscale video sample 
          sparams.width  = 176;   % width of the video
          sparams.height = 144;   % height of the video
          sparams.format = '420'; % YUV format ('420' for YUV 4:2:0 default)
        % sparams.sensmethod = 'binary'; % sensing matrix (multiple patterns) is binary
        sparams.sensmethod = 'mixhadamard'; % sensing matrix (multiple patterns) is binary
        sparams.SAMPLE_BINARY  = 0; % sample is binary (1-Y, 0-N)
        sparams.savesize = 5; % saving size of the image
        % sparams.binsize = 1; % binning size of the spatial light modulator
    REGENDATA = true; % flag of regenerating date (1-Y,0-N)
    % REGENDATA = false; % flag of regenerating date (1-Y,0-N)
    if REGENDATA
    % [1] generate and save dataset
        gendata3d(sparams);
    end

    % step 2. load and apply CS method to the dataset
    % [2] load dataset
    % [2.1] load simulation preferences
    load(sprintf('%s/sim_prefs.mat',datadir));
    % [2.2] load dataset
    load(sprintf('%s/%s%dby%d_samprate%.2f_snr%ddb.mat',simdatadir,...
        sparams.sampname,sparams.rows,sparams.nframe,sparams.samprate,...
        sparams.noisesnr));
    samp    = csdata.samp;    % sample
    sensmat = csdata.sensmat; % sensing matrixs
    meas    = csdata.meas;    % measurement vector

    rows    = sparams.rows;
    cols    = sparams.cols;
    nframe  = sparams.nframe;
    
    for imeth = 1:length(all3dmethod)
        % [3] apply CS method
        cs3dparams = []; % parameters for CS method\
            cs3dparams.rows     = rows;
            cs3dparams.cols     = cols;
            cs3dparams.nframe   = nframe;
            % cs3dparams.cs3dmethod = 'cs2d'; % 2D-CS for each frame
                % cs3dparams.srbasis  = 'haar'; % sparse representation basis
                % cs3dparams.srbasis  = 'daubechies'; % sparse representation basis
                cs3dparams.srbasis  = 'dct'; % sparse representation basis

                % cs3dparams.csmethod = 'gpsr'; % GPSR ell_1 solver
                % cs3dparams.csmethod = 'tv'; % TV regularization method
                cs3dparams.csmethod = 'gap'; % GAP ell_1 solver
                % cs3dparams.csmethod = 'gap-tv'; % GAP-TV  solver
                % cs3dparams.csmethod = 'gap2d'; % GAP2D solver
            % cs3dparams.cs3dmethod = 'tv3d'; % TV3D solver
            % cs3dparams.cs3dmethod = 'gap3d'; % GAP3D solver
            % cs3dparams.cs3dmethod = 'gap_tv3d'; % GAP_TV3D solver
            meth = all3dmethod{imeth};
            if strcmp(meth,'gap') || strcmp(meth,'tv')
                cs3dparams.cs3dmethod = 'cs2d'; % 2D-CS for each frame
                cs3dparams.csmethod = meth; % 2DCS method
            else
                cs3dparams.cs3dmethod = meth; % 2D-CS for each frame
            end

        profile clear
        profile -memory on
        tic;
        vout = cs3d(sensmat,meas,cs3dparams); % apply cs method
        t_cs = toc;
        profile report
        % profile -memory off

        % [4] save all the frames
        if strcmpi(cs3dparams.cs3dmethod,'cs2d')
            voutdir = sprintf('%s/%s%dby%dsamp%.1f/%s',outdir,sparams.sampname,rows,nframe,sparams.samprate,cs3dparams.csmethod);
        else
            voutdir = sprintf('%s/%s%dby%dsamp%.1f/%s',outdir,sparams.sampname,rows,nframe,sparams.samprate,cs3dparams.cs3dmethod);
        end
        if ~exist(voutdir,'dir')
            mkdir(voutdir);
        end
        for iframe = 1:nframe
            % sig_out = abs(sig_out);   % abs output
            % sig_out = max(sig_out,0); % threshold output
            samp_rc = reshape(vout(:,iframe),rows,cols);
            samp_nm = imnorm(samp_rc);
            samp_rz = imresize(samp_nm,sparams.savesize,'nearest');
            samp_or = reshape(samp(:,iframe),rows,cols); % original sample
            recrmse(iframe) = sqrt(immse(uint8(samp_nm*255),uint8(imnorm(samp_or)*255))); % root mean-square-error
            recpsnr(iframe) = psnr(uint8(samp_nm*255),uint8(imnorm(samp_or)*255));        % peak signal-to-noise-ratio
            recssim(iframe) = ssim(uint8(samp_nm*255),uint8(imnorm(samp_or)*255));        % structure similarity
            imwrite(samp_rz,sprintf('%s/%s%d_%02d_samprate%.2f_snr%ddb_%s_rmse%.2f_psnr%.2f_ssim%.2f_t%.1f.png',voutdir,...
                sparams.sampname,sparams.rows,iframe,sparams.samprate,sparams.noisesnr,cs3dparams.csmethod,recrmse(iframe),recpsnr(iframe),recssim(iframe),t_cs));
        end
    
    end % [loop] varying 3DCS methods
end % [loop] varying sampling rates

