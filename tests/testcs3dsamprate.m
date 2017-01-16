%% figure04 3DCS [PSNR] varying sampling rates
allsamprate = [10 20 40 60 80]*1e-2;
all3dmethod = {'tv','gap','tv3d','gap3d','gap_tv3d'};
allspec     = {'ko--','ks--','k*-','k+-','k^-'};
allpsnr = ...
[   16.91 	20.15 	24.36 	31.12 	38.02 ;
    13.18 	15.14 	17.86 	21.62 	11.78 ;
    19.60 	22.05 	30.09 	36.71 	42.61 ;
    18.47 	21.02 	18.98 	22.55 	22.79 ;
    12.33 	15.73 	27.95 	35.42 	42.59
];
allpsnrstd = ...
[   1.00 	1.74 	1.71 	1.39 	0.92 ;
    0.52 	0.72 	0.90 	0.61 	0.24 ;
    2.20 	3.25 	1.30 	1.87 	1.23 ;
    0.60 	0.68 	1.05 	0.96 	1.23 ;
    0.27 	0.93 	1.44 	1.60 	1.39
];
figure('position', [500, 500, 400, 300])  % create new figure with specified size  
m = length(all3dmethod);
for i = 1:m
    % errorbar(allsamprate,allpsnr(i,:),allpsnrstd(i,:),allspec{i}); hold on;
    plot(allsamprate,allpsnr(i,:),allspec{i}); hold on;
end
% ylim([0 56]);
xlabel('sampling rate');
ylabel('mean PSNR (dB)');
legend('TV (2D)','GAP (2D)','TV3D','GAP3D','GAP-TV (3D)','location','northwest');
legend boxoff;
% saveTightFigure('../../report/fig/fig05_psnr.pdf');
%% figure05 3DCS [SSIM] varying sampling rates
allsamprate = [10 20 40 60 80]*1e-2;
all3dmethod = {'tv','gap','tv3d','gap3d','gap_tv3d'};
allspec     = {'ko--','ks--','k*-','k+-','k^-'};
allssim = ...
[   0.52 	0.71 	0.88 	0.96 	0.99 ;
    0.18 	0.35 	0.60 	0.77 	0.22 ;
    0.70 	0.86 	0.96 	0.99 	1.00 ;
    0.62 	0.74 	0.81 	0.83 	0.84 ;
    0.22 	0.47 	0.93 	0.98 	0.99 
];
allssimstd = ...
[   0.03 	0.02 	0.01 	0.01 	0.00 ;
    0.02 	0.02 	0.01 	0.01 	0.01 ;
    0.03 	0.01 	0.00 	0.01 	0.01 ;
    0.02 	0.01 	0.01 	0.00 	0.01 ;
    0.01 	0.05 	0.01 	0.01 	0.01
];
figure('position', [500, 500, 400, 300])  % create new figure with specified size  
m = length(all3dmethod);
for i = 1:m
    % errorbar(allsamprate,allpsnr(i,:),allpsnrstd(i,:),allspec{i}); hold on;
    plot(allsamprate,allssim(i,:),allspec{i}); hold on;
end
% ylim([0 56]);
xlabel('sampling rate');
ylabel('mean SSIM');
legend('TV (2D)','GAP (2D)','TV3D','GAP3D','GAP-TV (3D)','location','northwest');
legend boxoff;
% saveTightFigure('../../report/fig/fig06_ssim.pdf');
%% figure07 3DCS [time complexity] varying sampling rates
allsamprate = [10 20 40 60 80]*1e-2;
all3dmethod = {'tv','gap','tv3d','gap3d','gap_tv3d'};
allspec     = {'ko--','ks--','k*-','k+-','k^-'};
alltime = ...
[  275	   283.1   318.2   404.5   605.8 ;
    10.1	19.5	27.6	57.2	79.5 ;
    76.1	84.6	65.8	78.2	78.9 ;
    11.4	22.1	28.7	42.8	83.8 ;
     5.9	11.6	21.2	31.8	42.5
];
% ylim([0 1]);
figure('position', [500, 500, 400, 300])  % create new figure with specified size  
m = length(all3dmethod);
for i = 1:m
    % errorbar(allsamprate,allpsnr(i,:),allpsnrstd(i,:),allspec{i}); hold on;
    semilogy(allsamprate,alltime(i,:),allspec{i}); hold on;
    % plot(allsamprate,alltime(i,:),allspec{i}); hold on;
end
% ylim([0 56]);
xlabel('sampling rate');
ylabel('running time (s)');
legend('TV (2D)','GAP (2D)','TV3D','GAP3D','GAP-TV (3D)','location','northwest');
legend boxoff;
% saveTightFigure('../../report/fig/fig07_time.pdf');