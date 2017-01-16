function [ vout ] = cs3d( sensmat,meas,params )
%CS3D Compressive sensing algorithms in solving three dimensional (3D)
%sparse reconstruction problem.
%   vout=CS3D(sensmat,meas,params) returns the N-by-F matrix of the sparse 
%   reconstruction result, where sensmat is the M-by-N-by-F sensing matrix 
%   or measurement matrix, M is the number of measurements, N is the 
%   Nyquist number of the orignal signal, F is the number of frames in the
%   volume, meas is the M-by-F measurement matrix that satisies the forward
%   model of three dimensional compressiver sampling (3D-CS), FM_3DCS, i.e.
%   for each frame meas(:,k)=sensmat(:,:,k)*vout(:,k)+nois(:,k), and nois 
%   here is the measurement noise on meas, params is a structure specifying
%   the parameters of the solver, including the algorithm (method) with 
%   the corresponding options.
%   See also FM_3DCS, CS, TV3D, GAP3D, GAP_TV3D.
% 
%   Copyright(C) 2016-2017 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>.
%   Last modified Jan 7, 2017.

% % [0] add *utils*, *algorithms* and *packages* to operation path, might 
% % be done out of this function block
% addpath('../utils')
% addpath('../algorithms')
% addpath(genpath('../packages'))

% [1] basic parameters
opts = []; % options for various algorithms

% [2] apply various compressive sensing algorithms
switch lower(params.cs3dmethod) % [params] csmethod
    case 'cs2d'    % [2.1] apply 2D-CS method to each frame of the whole 
                   %       volume
        % apply 2D-CS method to each frame
        for iframe = 1:params.nframe
            vout(:,iframe) = cs(sensmat(:,:,iframe),meas(:,iframe),params);
        end
        % [end] CS2D
        %
        %
    case 'tv3d'    % [2.2] three dimensional total variztion (3DTV) 
                   %       regularization for 3DCS
        % [2.2.1] generate TV3D gradient operator
        rows = params.rows;
        cols = params.cols;
        tv3dop = gentv3dop(cols,rows);
        % [2.2.2] options for TV3D
        opts = [];
            opts.mu0      = 1;
            opts.mumax    = 1e10;
            opts.maxiter  = 500;
            opts.rho      = 1.05;
            opts.tol      = 1e-4;
        lambda = 1;
        vout = tv3d(tv3dop,sensmat,meas,lambda,opts);
        % [end] TV3D
        % 
        % 
    case 'gap3d'   % [2.3] generalized alternating projection (GAP) in 
                   %       three dimensions (3D)
        % [2.3.1] expand the dimensions of rows and columns to the power of
        %         two (2^k)
        rows = params.rows;
        cols = params.cols;
        rk   = ceil(log2(rows));
        ck   = ceil(log2(cols));
        rows_exp = 2^rk; % expanded number of rows (2^rk)
        cols_exp = 2^ck; % expanded number of columns (2^ck)
        [M,N,F] = size(sensmat);
        if rows<rows_exp || cols<cols_exp
            for ifra = 1:F
                sensmat_exp(:,:,ifra) = [sensmat(:,:,ifra) zeros(M,rows_exp*cols_exp-N)];
            end
        else
            sensmat_exp = sensmat;
        end
        
        % [2.3.2] parameters configuration for GAP3D
        spbasis.space    = 'wavelet'; % transform for space, 'wavelet' or 'dct'
        spbasis.time     = 'dct';     % transform for spectrum, 'wavelet' or 'dct', dct is always used. I f we use wavelet, T need to be the power of 2. we can use no, means no transformation
        % spbasis.spectrum = 'no';  % Here we use no, means no transfromation in spectrum
        weighttype.space = 'tree';    % Here we can select:  'tree' or 'block'
        weighttype.time  = 'block';   % Here we can select:  'tree' or 'block', if we use tree, T need to be the power of 2
        weight_base.type = 'cos'; % here we can select 'exp' or 'cos'
        if strcmp(weight_base.type,'exp')
            weight_base.space = 1.5;   % This weight is the base of exponential decay. should be larger than 1 [1 2] is always used
            weight_base.time  = 1.5;
            weight_base.T     = 1.5;
        end
        % The block size for group
        block.row = 2;
        block.col = 2;
        block.T   = F/2;
        % stop criterion
        stopc.iternum = 50;
        stopc.err     = 10^-5;
        acc           = 2; % GAP with acceleration
        % acc           = 0; % GAP without acceleration
        ydim          = rows_exp*cols_exp;
        m_star = ceil(ydim/(block.row*block.col*block.T));
        
        A_3dcs  = @(x) fm_3dcs(sensmat_exp,x);
        At_3dcs = @(y) fmt_3dcs(sensmat_exp,y);
        
        % [2.3.3] apply GAP3D algorithm
        theta0 = gap3d(meas,A_3dcs,At_3dcs,rows_exp,cols_exp,F,block,...
                       spbasis,m_star,stopc,acc,weight_base,weighttype);
        vout = reshape(theta0(1:rows,1:cols,:),[rows*cols,F]);
        % [end] GAP3D
        % 
        % 
    case 'gap_tv3d' % [2.4] Total variation (TV)-based generalized 
                    %       alternating projection method in three 
                    %       dimensions (3D).
        % [2.4.1] options configuration for GAP_TV3D
        opts.rows     = params.rows;
        opts.cols     = params.cols;
        opts.nframe   = params.nframe;
        
        opts.lambda   = 1;
        opts.tvweight = 0.07;
        opts.tviter   = 5;
        opts.maxiter  = 50;
        opts.acc      = 1; % acceleration flag
        
        % [2.4.2] forward model of 3D-CS
        A_3dcs  = @(x) fm_3dcs(sensmat,x);
        At_3dcs = @(y) fmt_3dcs(sensmat,y);
        
        % [2.4.3] apply GAP_TV3D algorithm
        x_3d = gap_tv3d(meas,A_3dcs,At_3dcs,opts);
        vout = reshape(x_3d,[opts.rows*opts.cols,opts.nframe]);
        % [end] GAP_TV3D
        % 
        % 
    otherwise
        error('Unsupported three dimensional compressive sensing algorithm %s.\n',params.cs3dmethod);
end

