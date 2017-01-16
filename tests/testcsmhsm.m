%TESTCSMHSM test compressive sensing (CS) algorithms for sparse 
%reconstruction to form a table of the mixed Hadamard sensing matrix (MHSM)
%   See also CS.
% step 1. generate and save simulation dataset
clear; clc;
% close all;
allsamprate = [10 20 40 60 80]*1e-2;
allmethod   = {'tv','gap'};
allbasis    = {'dct'};

sparams = []; % paramters for simulation
for isamp = 1:length(allsamprate)
    sparams.samprate = allsamprate(isamp);
    % [1.0] parameters for simulation
        sparams.rows = 64; % number of rows (m)
        sparams.cols = 64; % number of columns (n)
        % sparams.samprate = 0.05; % sampling rate (gamma)
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

    % step 2. load and apply CS method to the dataset
    % [2] load dataset
    datadir    = '../../data';       % data directory
    simdatadir = '../../data/sim'; % simulation data directory
    outdir     = '../out'; % simulation output directory
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
    
    
    csparams = []; % parameters for CS method\
    for imeth = 1:length(allmethod)
        csparams.csmethod = allmethod{imeth};
        for ibas = 1:length(allbasis)
            csparams.srbasis = allbasis{ibas};
            % [3] apply CS method

                csparams.rows     = rows;
                csparams.cols     = cols;
                % csparams.srbasis  = 'haar'; % sparse representation basis
                % csparams.srbasis  = 'dct'; % sparse representation basis

                % csparams.csmethod = 'gpsr'; % GPSR ell_1 solver
                % csparams.csmethod = 'tv'; % TV regularization method
                % csparams.csmethod = 'gap'; % GAP ell_1 solver
                % csparams.csmethod = 'gap-tv'; % GAP ell_1 solver
            tic;
            sig_out = cs(sensmat,meas,csparams); % apply cs method
            t_cs = toc;
            samp_rc = reshape(sig_out,rows,cols);
            samp_nm = imnorm(samp_rc);
            samp_rz = imresize(samp_nm,sparams.savesize,'nearest');
            recrmse = sqrt(immse(uint8(samp_nm*255),uint8(imnorm(samp)*255))); % root mean-square-error
            recpsnr = psnr(uint8(samp_nm*255),uint8(imnorm(samp)*255));        % peak signal-to-noise-ratio
            recssim = ssim(uint8(samp_nm*255),uint8(imnorm(samp)*255));        % structure similarity
            if ~isfield(csparams,'srbasis')
                imwrite(samp_rz,sprintf('%s/%s%d_samprate%.2f_snr%ddb_%s_rmse%.2f_psnr%.2f_ssim%.2f_t%.1f.png',outdir,...
                    sparams.sampname,sparams.rows,sparams.samprate,sparams.noisesnr,csparams.csmethod,recrmse,recpsnr,recssim,t_cs));
            else
                imwrite(samp_rz,sprintf('%s/%s%d_samprate%.2f_snr%ddb_%s-%s_rmse%.2f_psnr%.2f_ssim%.2f_t%.1f.png',outdir,...
                    sparams.sampname,sparams.rows,sparams.samprate,sparams.noisesnr,csparams.csmethod,csparams.srbasis,recrmse,recpsnr,recssim,t_cs));
            end
            % figure; 
            % imshow(samp_rz);
        end
    end
end


%% figure [PSNR]
allpsnr = ...
[ 9.66 	12.68 	19.42 	20.83 	18.31;
  9.11	 9.92	13.60 	14.63 	16.11; 
 10.82 	18.40 	34.58 	40.71 	53.12;
  8.46	11.12	13.13 	17.10 	20.46 
];
figure;
plot(allsamprate,allpsnr(1,:),'ko--'); hold on;
plot(allsamprate,allpsnr(2,:),'ks--');
plot(allsamprate,allpsnr(3,:),'k*-');
plot(allsamprate,allpsnr(4,:),'k+-');
ylim([0 56]);
xlabel('sampling rate');
ylabel('PSNR');
legend('TV (binary)','GAP-DCT (binary)','TV (Hadamard)','GAP-DCT (Hadamard)');
legend boxoff;
% saveTightFigure('../../report/fig/fig06.pdf');
%% figure [SSIM]
allpsnr = ...
[   0.27 	0.47 	0.71 	0.69 	0.59;
    0.13 	0.19 	0.39 	0.47 	0.52;
    0.27 	0.69 	0.98 	0.99 	1.00;
    0.12	0.25	0.40 	0.55 	0.62 
];
figure;
plot(allsamprate,allpsnr(1,:),'ko--'); hold on;
plot(allsamprate,allpsnr(2,:),'ks--');
plot(allsamprate,allpsnr(3,:),'k*-');
plot(allsamprate,allpsnr(4,:),'k+-');
ylim([0 1]);
xlabel('sampling rate');
ylabel('SSIM');
legend('TV (binary)','GAP-DCT (binary)','TV (Hadamard)','GAP-DCT (Hadamard)');
legend boxoff;
% saveTightFigure('../../report/fig/fig07.pdf');

%% figure [time complexity]
allpsnr = ...
[   33.50 	33.90 	39.10 	50.70 	76.00;
     3.40 	 7.10 	13.40 	19.70 	26.20;
    34.40 	35.30 	41.10 	51.00 	75.20;
     3.50    7.10 	13.50 	20.70 	26.50 
];
figure;
plot(allsamprate,allpsnr(1,:),'ko--'); hold on;
plot(allsamprate,allpsnr(2,:),'ks--');
plot(allsamprate,allpsnr(3,:),'k*-');
plot(allsamprate,allpsnr(4,:),'k+-');
% ylim([0 1]);
xlabel('sampling rate');
ylabel('time (s)');
legend('TV (binary)','GAP-DCT (binary)','TV (Hadamard)','GAP-DCT (Hadamard)');
legend boxoff;
% saveTightFigure('../../report/fig/fig08.pdf');
