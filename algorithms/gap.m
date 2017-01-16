function [ xopt,errs ] = gap( A,y,opts )
%GAP generalized alternating projection algorithms in solving weighted
%ell_2,1 minimization problem
%   [xopt,errs]=GAP(A,y,opts) returns the optimal solution xopt and all the
%   errors errs during the iterations, where A is M-by-N sensing matrix, y
%   is M-by-1 measurement vector with Gaussian distributed noise, and opts
%   are options of the solver attached with some defaults.
%   [model] optimization model of the weighted ell_2,1 minimization problem
%     min  1/2*|w-theta|_2^2,
%     s.t. |theta|_{ell_{2,1}^{G,beta}}<=C,
%           A*w=y.
%   applying the alternating projection method of multipliers (ADMM)-like
%   method to the solver. Note that this method is different from
%   alternating projection in that the weighted ell_2,1 minimization method
%   is applied here rather than the ell_1 minimization method, and the
%   original ell_1 minimization problem can be written
%     min  |w|_1,
%     s.t.  A*x=y.
%   [default] parameters setting
%     opts.mstar                         $rows_of_sensing_matrix$
%     opts.tol                                               1e-4
%     opts.miniter                                             10
%     opts.maxiter                                            1e4
%   [references] This implementation follows the paper published on SIAM J. 
%   Imaging Sciences in 2014.
%     Liao, X., Li, H. & Carin, L. Generalized Alternating Projection for 
%       Weighted-$\ell_{2,1}$ Minimization with Applications to Model-Based
%       Compressive Sensing. SIAM Journal on Imaging Sciences 7, 797-823, 
%       doi:10.1137/130936658 (2014).
%   Copyright(C) 2016-2017 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>.
%   Last modified Dec 17, 2016.

% [0] default options
deft.mstar   = size(A,1); % [default] m_star
deft.tol     = 1e-4;      % [default] convergence tolerent
deft.miniter = 10;        % [default] minimum iterations
deft.maxiter = 1e3;       % [default] maximum iterations

opts = setdefault(opts,deft); % set default values to opts 

names = fieldnames(deft); % extract all the field names and assign them
for iname = 1:length(names)
    eval(sprintf('%s=opts.%s;',names{iname},names{iname}));
end

% [1] initialization
[~,N] = size(A);
theta = zeros(N,1);
u     = zeros(N,1);
yt    = y;
% [1.2] pre-calculation
pinvmat = A'/(A*A'); % right pseudo-inverse matrix

% [2] start iterations
conv = false; % convergence flag
for it = 1:maxiter
    % [2.1] update w_tilde^(t), y(t) and theta^(t)
    wt      = theta + pinvmat*(yt-A*theta);
    yt      = yt + (y-A*theta);
    lamball = sort(abs(wt),'descend');
    lambda  = lamball(mstar+1);
    theta   = sign(wt).*max(abs(wt)-lambda,0);
    % [2.2] update w^(t) and u^(t)
    w       = wt - u;
    u       = u + (w-theta);
    % [2.3] update errors
    errs(it) = norm(A*(theta-w),2);
    % [2.4] convergence criterion
    if it>=miniter && (max(abs(diff(errs(end-3:end)))/errs(end))<tol...
                        || max(errs(end-2:end))<tol)
        conv = true;
        break;
    end
end
if conv % converged at less than the maximum iterations
    fprintf('GAP solver converged at its %4d iterations.\n',it);
else
    fprintf('GAP solver reached the maximum %4d iterations.\n',it);
end    

% [3] output of the optimal solution and all the errors 
xopt = (w+theta)/2;

end

