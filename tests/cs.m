function [ sig_out ] = cs( sensmat,meas,params )
%CS compressive sensing algorithms in solving sparse reconstruction problem
%   out=CS(sensmat,meas,params) returns the N-by-1 vector of the sparse 
%   reconstruction result, where sensmat is the M-by-N sensing matrix or 
%   measurement matrix, M is the number of measurements and N is the 
%   Nyquist number of the orignal signal, meas is the M-by-1 measurement 
%   vector that satisies the compressiver sampling equation 
%   meas=sensmat*sig_out+nois, and nois here is the measurement noise on 
%   meas, params is a structure specifying the parameters of the solver, 
%   including the algorithm (method) with the corresponding options.
%   See also GPSR, TV, TVAL3, GAP_TV.
% 
%   Copyright(C) 2016-2017 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>.
%   Last modified Jan 7, 2017.

% % [0] add *algorithms* and *packages* to operation path, might be done out
% % of this function block
% addpath('../algorithms')
% addpath(genpath('../packages'))

% [1] basic parameters
opts = []; % options for various algorithms

% [2] apply various compressive sensing algorithms
switch lower(params.csmethod) % [params] csmethod
    case 'gpsr'     % [2.1] gradient projection sparse reconstruction
        % [2.1.1] sparse representation
        rows = params.rows;
        cols = params.cols;
        switch lower(params.srbasis) % [params] srbasis
            case 'haar'    % Haar wavelet transform basis
                Fr = haarmtx(rows); % assist rows is power of 2
                if cols==rows
                    Fc = Fr;
                else
                    Fc = haarmtx(cols); % assist cols is power of 2
                end
                Psi = kron(Fc,Fr);  % sparsifing basis
            case 'dct'     % discrete cosine transform basis
                Dr = dctmtx(rows); % assist rows is power of 2
                if cols==rows
                    Dc = Dr;
                else
                    Dc = dctmtx(cols); % assist cols is power of 2
                end
                Psi = kron(Dc,Dr);  % sparsifing basis
            otherwise
                error('unsupported sparse representation basis %s.\n',params.srbasis);
        end
        A = double(sensmat)*Psi;
        
        % [2.1.2] ell_1 solver gpsr 
        % parameters configuration of gpsr
%         [M,N] = size(sensmat);
        sig0 = zeros(size(sensmat,2),1); % start point
        lambda = 1; % regularization parameter for ell_1-ell_2 optimization
%         sig0 = params.sig0; % start point
%         lambda = params.lambda; % regularization parameter for ell_1-ell_2 optimization
        
        tol = 1e-4; % convergence tolerance for gradient projection method
        opts.mu      = 0.1;   % [backtracking l.s.] (0,1/2)
        opts.beta    = 0.5;   % [backtracking l.s.] (0,1)
        opts.maxiter = 1000; % [iteration] maximum iteration
        % apply GPSR algorithm
        w = gpsr(A,meas,sig0,lambda,tol,opts);
        sig_out = Psi*w;
        % [end] GPSR
        %
        % 
    case 'tv'      % [2.2] total variation regularization (different from tval3)
        % [2.2.0] parameters configuration of tv
        rows = params.rows;
        cols = params.cols;
        % [2.2.1] generator tv operator (gradient operator)
        gradord = 1; % gradient order
        tvop = gentvop(rows,cols,gradord); 
        sensmat = double(sensmat);
        % [2.2.2] apply tv solver
        opts.x0      = pinv(sensmat)*meas; % start point (pseudo-inverse result)
        opts.mu0     = 1;    % initial value of penalty factor mu
        opts.mubar   = 1e10; % maximum value of penalty factor mu
        opts.rho     = 1.05; % multiplication step of penalty factor mu
        opts.tol     = 1e-4;
        opts.miniter = 20;   % minimum iterations
        opts.maxiter = 500;  % maximum iterations
        % apply TV algorithm
        sig_out = tv(sensmat,meas,tvop,opts);
        % [end] TV
        % 
        % 
    case 'gap'     % [2.3] generalized alternating projection (GAP)
                   %       applying mixed Hadamard sensing method (MHSM)
        rows = params.rows;
        cols = params.cols;
        switch lower(params.srbasis) % [params] srbasis
            case 'haar'       % Haar wavelet transform basis
                Fr = haarmtx(rows); % assist rows is power of 2
                if cols==rows
                    Fc = Fr;
                else
                    Fc = haarmtx(cols); % assist cols is power of 2
                end
                Psi = kron(Fc,Fr);  % sparsifing basis
            case 'daubechies' % Daubechies wavelet transform basis (DB-8)
                level = 3;
                qmf = MakeONFilter('Daubechies',8); % Daubechies-8 wavelet
                sig_level_row = log2(rows);
                sig_level_col = log2(cols);
                Fr = get_waveletMatrix(qmf,sig_level_row,level,level);
                Fc = get_waveletMatrix(qmf,sig_level_col,level,level);
                Psi = kron(Fc,Fr);  % sparsifing basis
            case 'dct'        % discrete cosine transform basis
                Dr = dctmtx(rows); % assist rows is power of 2
                if cols==rows
                    Dc = Dr;
                else
                    Dc = dctmtx(cols); % assist cols is power of 2
                end
                Psi = kron(Dc,Dr);  % sparsifing basis
            otherwise
                error('unsupported sparse representation basis %s.\n',params.srbasis);
        end
        A = double(sensmat)*Psi;
        
        opts.tol     = 1e-3;
        opts.maxiter = 300;
        % apply GAP algorithm
        w = gap(A,meas,opts);
        sig_out = Psi*w;
        % [end] GAP
        % 
        % 
    case 'gap-tv'  % [2.4] generalized alternating projection (GAP) in 
                   % solving total variation (TV) minimization problem
        % [2.4.0] parameters configuration of tv
        rows  = params.rows;
        cols  = params.cols;
        % [2.4.1] generator tv operator (gradient operator)
        gradord = 1; % gradient order
        tvop = gentvop(rows,cols,gradord); 
        sensmat = double(sensmat);
        
        % [2.4.2] apply GAP-TV solver
        opts.lambda  = 16;        % TV regulizer
        opts.alpha   = 8;         % parameter for clipping algorithm, 
                                  %   alpha>=max(eig(tvop*tvop')) [=8]
        opts.tol     = 1e-5;      % convergence tolerent
        opts.acc     = false;     % accelerated version of GAP method
        % apply GAP_TV algorithm 
        sig_out = gap_tv(sensmat,meas,tvop,opts);
        % [end] GAP-TV
        % 
        % 
    case 'gap2d'   % [2.5] generalized alternating projection (GAP) for
                   %       two-dimensional (2D) compressive sensing
        % [2.5.1] expand the dimensions of rows and columns to be the power
        %         of two (2^k)
        rows = params.rows;
        cols = params.cols;
        rk   = ceil(log2(rows));
        ck   = ceil(log2(cols));
        rows_exp = 2^rk; % expanded number of rows (2^rk)
        cols_exp = 2^ck; % expanded number of columns (2^ck)
        [M,N] = size(sensmat);
        semsmat_exp = [sensmat zeros(M,rows_exp*cols_exp-N)];
        
        % [2.5.2] parameters configuration of GAP2D
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
        % block.T   = T/2;
        block.T   = 1; % for 2D-CS
        % stop criterion
        stopc.iternum = 50;
        stopc.err     = 10^-5;
        acc           = 2; % GAP with acceleration
        % acc           = 0; % GAP without acceleration
        ydim          = rows_exp*cols_exp;
        m_star = ceil(ydim/(block.row*block.col*block.T));
        
        A_cs  = @(x) fm_cs(semsmat_exp,x);
        At_cs = @(y) fmt_cs(semsmat_exp,y);
        
        T = 1; % number of frames varying time
        m_star = m_star-1; % m_star should be decreased for 2D-CS
        % [2.5.3] apply GAP2D algorithm
        theta0 = gap2d(meas,A_cs,At_cs,rows_exp,cols_exp,T,block,...
                       spbasis,m_star,stopc,acc,weight_base,weighttype);
        
        theta_vec = theta0(:);
        sig_out = theta_vec(1:rows*cols); % recover the reconstructed signal
        % [end] GAP2D
        % 
        % 
    case 'gap3d'   % [2.6] generalized alternating projection (GAP) in 
                   %       three dimensions (3D)
        % 
        % to apear in CS3D
        % 
        % 
        % [end] GAP3D
        % 
        % 
    otherwise
        error('Unsupported compressive sensing algorithm %s.\n',params.csmethod);
end

