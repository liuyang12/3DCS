function [ x_opt ] = gpsr( A,y,x0,lambda,tol,opts )
%GRSP gradient projection for sparse reconstruction
%   x_opt=GPSR(A,y,tol) returns the optimal sparse solution for ell_1-ell_2
%   optimization appling gradient projection method, where A is the m-by-n 
%   matrix, y is the m-by-1 vector, x0 is the start point,lambda is the 
%   regularization parameter for ell_1-ell_2 optimization, tol is the 
%   convergence tolerance for gradient projection method, and opts are 
%   the options for the algorithm.
%   the ell-1_ell_2 optimization problem goes to be:
%   min_x 1/2*||y-A*x||_{ell_2}^2+lambda*||x||_{ell_1}.
mu      = opts.mu;      % [backtracking l.s.] (0,1/2)
beta    = opts.beta;    % [backtracking l.s.] (0,1)
maxiter = opts.maxiter; % [iteration] maximum iteration
alpha_min = 0.01;
alpha_max = 1;
[m,n] = size(A);
% transfer the ell_1-ell-2 optimization problem to quadratic problem with
% inequality constrained.
%   min_z f(z)=1/2*z'*B*z+c'*z  s.t. z>=0.
%   nabla f(z)=B*z+c.
u0 = max(0,x0);
v0 = u0-x0;
z0 = [u0;v0];
b = A'*y;
c = lambda+[-b;b];
X = A'*A;
B = [X,-X;-X,X];

% start iteration
zk = z0;
conv = false;
for k = 0:maxiter
    Dfzk = B*zk+c; % gradient at zk
    gk = zeros(2*n,1);
    pospointind = find(zk>0); % positive point index
    neggradind = find(Dfzk(pospointind)<0); % negative gradient index 
    gk(pospointind(neggradind)) = Dfzk(pospointind(neggradind));
    alpha = gk'*gk/(gk'*B*gk);
    alpha = max(alpha_min,alpha);
    alpha = min(alpha_max,alpha);
    while 1
        zpk = max(0,zk-alpha*Dfzk);
        if (1/2*zpk'*B*zpk+c'*zpk<=1/2*zk'*B*zk+c'*zk-mu*Dfzk'*(zk-zpk))
            conv = true;
            break;
        else
            alpha=alpha*beta;
        end
    end
    zk0 = zk;
    zk = zpk;
    if norm(min(zk,B*zk+c))<=tol || norm(zk-zk0)<=tol
        conv = true;
        break;
    end
end
if conv % converged at less than the maximum iterations
    fprintf('GPSR solver converged at its %4d iterations.\n',k);
else
    fprintf('GPSR solver reached the maximum %4d iterations.\n',k);
end 

x_opt = zk(1:n)-zk(n+1:end);

end

