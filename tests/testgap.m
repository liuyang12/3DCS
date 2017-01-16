%TESTGAP test generalized alternating projection (GAP) solver for ell_1
%minimization problem
%   See also GAP.
clear; clc;
% close all;
% simulating data
allsamprate = [.5 1 2 5 10 20 40]*1e-2; % all the sampling rates

reptime  = 100;    % repeat sampling times
N        = 1000;   % Nyquist number
samprate = 0.5e-2; % sampling rate
sparsity = 40;     % sparsity
noisesnr = 50;     % SNR (dB)

M = round(samprate*N);
k = sparsity;
A = randn(M,N); 
x0 = zeros(N,1);
idx = randperm(N);
i1 = idx(1:k); 
x0(i1) = randn(k,1);
y0 = A*x0;
noise = randn(M,1);
noise = noise*10^(-noisesnr/20)*sqrt(y0'*y0)/sqrt(noise'*noise);
snr = 10*(log10(y0'*y0) - log10(noise'*noise));
y  = y0 + 0*noise;     
    
% applying GAP to reconstruct the sparse signal
opts = [];
opts.tol = 1e-4;
tic; 
[xopt_gap,errs] = gap(A,y,opts);  
t_gap = toc;
mse_gap = immse(xopt_gap,x0);
sc_gap = sparcons(xopt_gap,x0,opts.tol/10);

% applying GPSR to reconstuct the sparse signal
sig0 = zeros(N,1); % start point
lambda = 1; % regularization parameter for ell_1-ell_2 optimization
tol = 1e-4; % convergence tolerance for gradient projection method
opts.mu      = 0.1;   % [backtracking l.s.] (0,1/2)
opts.beta    = 0.5;   % [backtracking l.s.] (0,1)
opts.maxiter = 1000; % [iteration] maximum iteration
tic;
xopt_gpsr = gpsr(A,y,sig0,lambda,tol,opts);
t_gpsr = toc;
mse_gpsr = immse(xopt_gpsr,x0);
sc_gpsr = sparcons(xopt_gpsr,x0,tol/10);

% plotting results
figure;
subplot(3,1,1); 
stem(x0,'k+'); ylim([-2 2]);
title('original signal');
subplot(3,1,2); 
stem(xopt_gap,'k+'); ylim([-2 2]);
title(sprintf('GAP method (MSE = %.1e, SC = %.1f, t = %.2f s)',mse_gap,sc_gap,t_gap));
subplot(3,1,3); 
stem(xopt_gpsr,'k+'); ylim([-2 2]);
title(sprintf('GPSR method (MSE = %.1e, SC = %.1f, t = %.2f s)',mse_gpsr,sc_gpsr,t_gpsr));
% saveTightFigure('../../report/fig/fig01.pdf');

figure; semilogy(errs,'-o'); grid; 
title('convergence of reconstruction error on linenar manifold'); 
ylabel('squared reconstruction error');
xlabel('iterations');
%figure(3); semilogy(model.g); grid; 
%title(sprintf('fitting error, using noisy measurements with SNR=%0.2f dB',snr),'FontSize',15);
% set(gca,'FontSize',20);