function [ xopt,errs ] = gap_tv( A,y,tvop,opts )
%GAP_TV generalized alternating projection (GAP) method in solving total
%variation (TV) minimizarion problem
%   [xopt,errs]=GAP_TV(A,y,tvop,opts) returns the optimal solution for 
%   GAP-TV solver and the errors in each iteration, where A is M-by-N 
%   sensing matrix, y is M-by-1 measurement vector with Gaussian 
%   distributed boise, tvop is the TV operation matrix and opts are options
%   of the solver attached with some defaults.
%   [model] optimization model of the total variation (TV) minimization
%   problem
%     min  1/2*|x-theta|_2^2+lambda*|TV(theta)|,
%     s.t. A*x=y.
%   applying the generalized alternating projection (GAP) method. Note that
%   GAP is for weighted ell_2,1 minimization purpose and each element of 
%   TV, TV(x) is regarded as a group here.
%   [default] parameters setting
%     opts.x0                       $pseudo-inverse_result$
%     opts.lambda                                         1
%     opts.alpha                                          8
%     opts.tol                                         1e-4
%     opts.acc                                         true
%     opts.miniter                                       10
%     opts.maxiter                                      1e4  
%   [refrences] This implementation follows the paper pubulished on Image 
%   Processing (ICIP), 2016 IEEE International Conference on.
%     Yuan, X. Generalized alternating projection based total variation 
%       minimization for compressive sensing, in Image Processing (ICIP), 
%       2016 IEEE International Conference on.  2539-2543.
%   Copyright(C) 2016-2017 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>.
%   Last modified Dec 17, 2016.

% [0] default options
deft.x0      = pinv(A)*y; % [default] start point (pseudo-inverse result)
deft.lambda  = 1;         % [default] regulizer
deft.alpha   = 8;         % [default] parameter for clipping algorithm, 
                          %   alpha>=max(eig(tvop*tvop'))
deft.tol     = 1e-4;      % [default] convergence tolerent
deft.acc     = true;      % [default] accelerated version of GAP method
deft.miniter = 10;        % [default] minimum iterations
deft.maxiter = 1e4;       % [default] maximum iterations

opts = setdefault(opts,deft); % set default values to opts 

names = fieldnames(deft); % extract all the field names and assign them
for iname = 1:length(names)
    eval(sprintf('%s=opts.%s;',names{iname},names{iname}));
end      

alpha_ = opts.alpha;

% [1] initialization
[~,N] = size(A);
% theta = zeros(N,1);
theta = x0;
z     = zeros(size(tvop,1),1);
yt    = y;
% [1.2] pre-calculation
pinvmat = A'/(A*A'); % right pseudo-inverse matrix

% [2] iterations
conv = false; % convergence flag
for it = 1:maxiter
    % [2.1] update x^(t), z^(t) and theta^(t)
    x     = theta + pinvmat*(yt-A*theta);
    if opts.acc % accelerated GAP
        yt = yt + (y-A*theta);
    end
    theta = x - tvop'*z;
    z     = clip(z+1/alpha_*tvop*theta,lambda/2);
    % [2.2] update errors
%     errs(it) = norm(tvop*(x+theta)/2,1);
    errs(it) = norm(A*(theta-x),2);
    % [2.3] convergence criterion
    if it>=miniter && (max(abs(diff(errs(end-3:end)))/errs(end))<tol...
                        || max(errs(end-2:end))<tol)
        conv = true;
        break;
    end
end
if conv % converged at less than the maximum iterations
    fprintf('GAP-TV solver converged at its %4d iterations.\n',it);
else
    fprintf('GAP-TV solver reached the maximum %4d iterations.\n',it);
end    

% [3] output of the optimal solution and all the errors 
xopt = (x+theta)/2;
end

% clipping function with thresold of T
function z = clip(b,T)
z = max(min(b,T),-T);
end

