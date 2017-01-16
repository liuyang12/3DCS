%TESTELL1SPARSITY test ell_1 minimization problem with respect to sparsity
%   See also GAP, GPSR.
clear; clc;
% close all;
% simulating data
allsparsity = [5 10 20 40 60 80 100]*1e-2; % all the sampling rates

reptime  = 5;    % repeat sampling times
N        = 1000;   % Nyquist number
samprate = 20e-2; % sampling rate
sparsity = 40;     % sparsity
noisesnr = 50;     % SNR (dB)

for ispar = 1:length(allsparsity)
    sparsity = allsparsity(ispar);
    for irep = 1:reptime
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
        t_gap(irep,ispar) = toc;
        mse_gap(irep,ispar) = immse(xopt_gap,x0);
        sc_gap(irep,ispar) = sparcons(xopt_gap,x0,opts.tol/10);

        % applying GPSR to reconstuct the sparse signal
        sig0 = zeros(N,1); % start point
        lambda = 1; % regularization parameter for ell_1-ell_2 optimization
        tol = 1e-4; % convergence tolerance for gradient projection method
        opts.mu      = 0.1;   % [backtracking l.s.] (0,1/2)
        opts.beta    = 0.5;   % [backtracking l.s.] (0,1)
        opts.maxiter = 1000; % [iteration] maximum iteration
        tic;
        xopt_gpsr = gpsr(A,y,sig0,lambda,tol,opts);
        t_gpsr(irep,ispar) = toc;
        mse_gpsr(irep,ispar) = immse(xopt_gpsr,x0);
        sc_gpsr(irep,ispar) = sparcons(xopt_gpsr,x0,tol/10);
    end
end

%% Sparse signal recovery results concerning accuracy varying sampling rates.
figure; subplot(121);
errorbar(allsparsity*100,mean(log10(mse_gap)),std(log10(mse_gap),1),'k+-'); hold on;
errorbar(allsparsity*100,mean(log10(mse_gpsr)),std(log10(mse_gpsr),1),'ko--');
title('mean squared error (MSE)');
xlabel('sampling rate (%)');
ylabel('log(MSE)');
legend('GAP','GPSR'); legend boxoff;

subplot(122);
errorbar(allsparsity*100,mean(sc_gap),std(sc_gap,1),'k+-'); hold on;
errorbar(allsparsity*100,mean(sc_gpsr),std(sc_gpsr),'ko--');
title('sparse consistancy (SC)')
xlabel('sampling rate (%)');
ylabel('SC');
legend('GAP','GPSR'); legend boxoff;
% saveTightFigure('../../report/fig/fig02.pdf');

%% Sparse signal recovery results concerning time complexity varying sampling rates.
figure;
errorbar(allsparsity*100,mean(t_gap),std(t_gap,1),'k+-'); hold on;
errorbar(allsparsity*100,mean(t_gpsr),std(t_gpsr),'ko--');
title('time complexity');
xlabel('sampling rate (%)');
ylabel('time (s)');
legend('GAP','GPSR'); legend boxoff;
% saveTightFigure('../../report/fig/fig03.pdf');
