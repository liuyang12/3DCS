%TESTCS3DEXP Test compressive sensing (CS) algorithms for three dimensional
%sparse reconstruction of experimental datasets.
%   See also CS3D, TESTCS3D.

%% load and apply CS method to the experimental dataset
% [1] load dataset
clear; clc;
% close all;
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

datadir    = '../../data';                   % data directory
simdatadir = '../../data/exp';               % experimental data directory
gtruthdir  = '../../data/exp/groundtruth';   % experimental ground truth
outdir     = '../expout';                    % experimental output directory
if ~exist(outdir,'dir')
    mkdir(outdir);
end

sparams.sampname = 'rman'; % [exp] sample name

% [1.1] load dataset
load(sprintf('%s/%s.mat',simdatadir,sparams.sampname));
% sensmat = P_8ch;    % sensing matrixs
% meas    = Y_8ch;    % measurement vector

[M,N,F] = size(sensmat);
rows    = round(sqrt(N));
cols    = rows;
nframe  = F;
sparams.samprate = M/N;
sparams.savesize = 4;

% [2] apply CS method
cs3dparams = []; % parameters for CS method
    cs3dparams.rows     = rows;
    cs3dparams.cols     = cols;
    cs3dparams.nframe   = nframe;
    cs3dparams.cs3dmethod = 'cs2d'; % 2D-CS for each frame
        % cs3dparams.srbasis  = 'haar'; % sparse representation basis
        % cs3dparams.srbasis  = 'daubechies'; % sparse representation basis
        cs3dparams.srbasis  = 'dct'; % sparse representation basis

        % cs3dparams.csmethod = 'gpsr'; % GPSR ell_1 solver
        cs3dparams.csmethod = 'tv'; % TV regularization method
        % cs3dparams.csmethod = 'gap'; % GAP ell_1 solver
        % cs3dparams.csmethod = 'gap-tv'; % GAP-TV  solver
        % cs3dparams.csmethod = 'gap2d'; % GAP2D solver
    % cs3dparams.cs3dmethod = 'tv3d'; % TV3D solver
    % cs3dparams.cs3dmethod = 'gap3d'; % GAP3D solver
    % cs3dparams.cs3dmethod = 'gap_tv3d'; % GAP_TV3D solver

profile clear
profile -memory on
tic;
vout = cs3d(double(sensmat),meas,cs3dparams); % apply cs method
t_cs = toc;
profile report
% profile -memory off

% [3] save all the frames
if strcmpi(cs3dparams.cs3dmethod,'cs2d')
    voutdir = sprintf('%s/%s%dby%d/%s',outdir,sparams.sampname,rows,nframe,cs3dparams.csmethod);
else
    voutdir = sprintf('%s/%s%dby%d/%s',outdir,sparams.sampname,rows,nframe,cs3dparams.cs3dmethod);
end
if ~exist(voutdir,'dir')
    mkdir(voutdir);
end
% % [3.1] load the groud truth of the dataset
% filenames = getallfiles(gtruthdir);
% for ifile = 1:length(filenames)
%     gtfile = filenames{ifile};
%     gtimg  = imread(gtfile);
%     samp(:,:,ifile) = imresize(gtimg,[rows cols],'nearest');
% end
    
for iframe = 1:nframe
    % vout = abs(vout);   % abs output
    vout = max(vout,0); % threshold output
    samp_rc = reshape(vout(:,iframe),rows,cols);
    samp_nm = imnorm(samp_rc);
    samp_rz = imresize(samp_nm,sparams.savesize,'nearest');
    % samp_or = reshape(samp(:,iframe),rows,cols); % original sample
    samp_or = imnorm(double(samp(:,:,iframe))); % original sample
    recrmse(iframe) = sqrt(immse(uint8(samp_nm*255),uint8(samp_or*255))); % root mean-square-error
    recpsnr(iframe) = psnr(uint8(samp_nm*255),uint8(samp_or*255));        % peak signal-to-noise-ratio
    recssim(iframe) = ssim(uint8(samp_nm*255),uint8(samp_or*255));        % structure similarity
    imwrite(samp_rz,sprintf('%s/%s%d_%02d_samprate%.2f_%s_rmse%.2f_psnr%.2f_ssim%.2f_t%.1f.png',voutdir,...
        sparams.sampname,rows,iframe,sparams.samprate,cs3dparams.csmethod,recrmse(iframe),recpsnr(iframe),recssim(iframe),t_cs));
end
figure; 
imshow(samp_rz);


