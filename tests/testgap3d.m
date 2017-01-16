%TESTGAP3D Test three dimensional (3D) generalized alternating projection 
%(GAP) algorithm, that is GAP3D algorithm.
clear; clc;
% close all;
% [0] add dependencies to path
addpath('../utils')
addpath('../algorithms')
addpath(genpath('../packages'))

% [1] load the CASSI dataset
datadir    = '../../data';       % data directory
cassidir   = '../../data/cassi'; % cassi data directory
cassipath  = sprintf('%s/CASSI_bird.mat',cassidir); % cassi data path
load(cassipath);

Phi = squeeze(Phi(:,:,:,1)); % Nx*Ny*Nt patterns
y   = Y(:,:,1); % Nx*Ny coded image
x0  = slit; % Nx*Ny*Nt original image sequence (ground truth),
            % that is the forward model y=sum(Phi.*x0,3)
[Nx,Ny,Nt] = size(Phi);

%% [2] obtain the linear forward model of the reconstruction problem
A_fun  = @(f) fm_cassi(Phi,f);
At_fun = @(y) fmt_cassi(Phi,y);

x_rec_gap3d = run_gap_cassi(y,Phi);

