function [ xopt ] = tv( A,y,H,opts )
%TV Total variation (TV) regularization for sparse reconstruction problem
%   xopt=TV(A,y,H,opts) returns the optimal TV regularization result, where
%   A is m-by-n measurement matrix, y is m-by-1 measurements with Gaussian
%   noise, H is the *-by-n TV operator, and opts are the parameters used in
%   this solver.
%   [model] optimization model of the sparse reconstruction solver
%     min  |Hx|_1,
%     s.t.  A*x=y.
%   by solving the augmented Lagrangian applying alternating direction
%     L(w,x)=|w|_1+mu/2*|H*x-w|_2^2+mu/2*|A*x-y|_2^2-r'*(H*x-w)-s'*(A*x-y).
%   note that the augmented Lagrangian is different from TV/L2 in TVAL3 by
%   puting the same penalty factor for |H*x-w|_2 and |A*x-y|_2.
%   [default] parameters setting:
%     opts.x0                     $pseudo-inverse_result$
%     opts.mu0                                          1
%     opts.mubar                                     1e10
%     opts.rho                                       1.05
%     opts.tol                                       1e-5
%     opts.miniter                                     20
%     opts.maxiter                                    500
%   [references] This implementation follows the paper published on SIAM J.
%   Imaging Sciences in 2008.
%     Wang, Y., Yang, J., Yin, W. & Zhang, Y. A New Alternating 
%       Minimization Algorithm for Total Variation Image Reconstruction. 
%       SIAM Journal on Imaging Sciences 1, 248-272, doi:10.1137/080724265 
%       (2008).
%   Copyright(C) 2016-2017 by <a href="matlab: 
%   web('https://github.com/liuyang12')">Y. Liu</a>.
%   Last modified Dec 17, 2016.

% [0] default options
deft.x0      = pinv(A)*y; % [default] start point (pseudo-inverse result)
deft.mu0     = 1;         % [default] initial value of penalty factor mu
deft.mubar   = 1e10;      % [default] maximum value of penalty factor mu
deft.rho     = 1.05;      % [default] multiplication step of penalty factor mu
deft.tol     = 1e-4;      % [default] convergence tolerance
deft.miniter = 20;        % [default] minimum iterations
deft.maxiter = 500;       % [default] maximum iterations

opts = setdefault(opts,deft); % set default values to opts

names = fieldnames(deft); % extract all the field names and assign them
for iname = 1:length(names)
    eval(sprintf('%s=opts.%s;',names{iname},names{iname}));
end

% [1.0] initialization
x = x0; % x0 start point
w = H*x; % sparse vector in gradient 
r = zeros(size(w)); % r Lagrange multiplier
s = zeros(size(y)); % s Lagrange multiplier
obj_k = norm(w(:),1); % object function
mu = mu0; % mu penalty factor

% [1.1] pre-calculation
pinvmat = pinv(H'*H+A'*A); % pseudo-inverse matrix

% [2] start iteration
conv = false; % convergence flag
for it = 0:maxiter % number of iteration count
    % [2.1] update optimization variables w and x
    intv = H*x-1/mu*r;
    fpos = intv> 1/mu;
    fneg = intv<-1/mu;
    w = (intv-1/mu).*fpos+(intv+1/mu).*fneg;
    x = pinvmat*(H'*(w+1/mu*r)+A'*(y+1/mu*s)); % sign of y-1/mn*s is wrong
    
    % [2.2] update Lagrange multipliers r and s
    r = r-mu*(H*x-w);
    s = s-mu*(A*x-y);
     
    % [2.3] update the penalty factor mu
    mu = min(mu*rho,mubar);
    
    % [2.4] convergence criterion
    obj_k0 = obj_k; % previous value of object function
    obj_k = norm(H*x,1); % updated value of object function
    crit1 = obj_k0-obj_k; % criterion #1
    crit2 = norm(obj_k0-obj_k)/norm(obj_k); % criterion #2
    if ((crit1<tol&&crit1>=0) || abs(crit1)<tol/100)...
            && crit2<tol... % object function 
            && it>miniter   % minimum iterations
        conv = true;
        break;
    end    
end
if conv % converged at less than the maximum iterations
    fprintf('TV solver converged at its %4d iterations.\n',it);
else
    fprintf('TV solver reached the maximum %4d iterations.\n',it);
end

% [3] output of the optimal solution
xopt = x; % optimal solution

end

