function X = tv3d(H,P,Z,lambda,opts)
%TV3D  Three dimensional total variation (3DTV) regularization for three
%dimensional (3D) compressive sensing reconstruction problem.
%   [model]
%   ---------------------------------------
%   min \sum_i |H*X_i|_1 + lambda*|X*OP|_1
%   s.t. P_i*X_i=Z_i, i = 1,...,F
%   ------------------------------------------
%   [variables]
%     N    -- Nyquist number
%     F    -- frame number
%     M    -- measurement number
%     P    -- measurement matrix array: [M, N, F]
%     Z    -- meassurement array: [M, F]
%     X_0  -- initial X: [N, F]
%     X_gt -- groundtrue X: [N, F]
%     H    -- gradient operator
%   See also GENTV3DOP, TV, CS3D.

mu0     = opts.mu0;
mumax   = opts.mumax;
maxiter = opts.maxiter;
rho     = opts.rho;
tol     = opts.tol;

[~,N,F]  = size(P);

%%% OP -- intra-frame operator
OP = -triu(ones(F,F-1),0)+triu(ones(F,F-1),-1)...
            -triu(ones(F,F-1),0)+triu(ones(F,F-1),1);

%%% lambda
% lambda=norm(H*double(X_0),1)/norm(double(X_0)*OP,1); 

%%% initialization
% X=X_0;
X = zeros(N, F);
mu = mu0;
W = H*X;
V = X*OP;
R = zeros(size(W));
S = zeros(size(V));
T = zeros(size(Z));
    
%%% 
% disp('precalculation of TV3D')
HT  = H';
HTH = H'*H;
I   = eye(N);
X_update_pinv = zeros(N,N,F);
for i = 1:F
    PTP = P(:,:,i)'*P(:,:,i);
    if i==1 || i==F
        X_update_pinv(:,:,i) = (HTH+I+PTP)\I;
    else
        X_update_pinv(:,:,i) = (HTH+2*I+PTP)\I;
    end
end
clear HTH PTP I

%% iteration
iter = 0;
conv = false; % convergence flag
object_kp = norm(H*X,1)+lambda*norm(X*OP,1);
% disp('main iterations of TV3D');
while ~conv
    iter = iter+1;
    % update X
    for i = 1:F
        if i==1
            X(:,i) = X_update_pinv(:,:,i)...
                        *(  HT*(W(:,i)+1/mu*R(:,i))...
                        +(  X(:,i+1)-V(:,i)-1/mu*S(:,i)  )...
                        + P(:,:,i)'*(Z(:,i)-1/mu*T(:,i)) );
        elseif i==F
            X(:,i) = X_update_pinv(:,:,i)...
                        *(  HT*(W(:,i)+1/mu*R(:,i))...
                        +(  X(:,i-1)+V(:,i-1)+1/mu*S(:,i-1)  )...
                        + P(:,:,i)'*(Z(:,i)-1/mu*T(:,i)) );            
        else
            X(:,i) = X_update_pinv(:,:,i)...
                        *(  HT*(W(:,i)+1/mu*R(:,i))...
                        +( X(:,i+1)-V(:,i)-1/mu*S(:,i)...
                        +  X(:,i-1)+V(:,i-1)+1/mu*S(:,i-1)  )...
                        + P(:,:,i)'*(Z(:,i)-1/mu*T(:,i)) );              
        end
    end  

    % update W
    tmp = H*X-1/mu*R;
    flag1 = (tmp>1/mu);
    flag2 = (tmp<-1/mu);
    W = (tmp-1/mu).*flag1+(tmp+1/mu).*flag2;  

    % update V
    tmpV = X*OP-S/mu;
    flag1 = (tmpV>lambda/mu);
    flag2 = (tmpV<-lambda/mu);
    V = (tmpV-lambda/mu).*flag1+(tmpV+lambda/mu).*flag2;  % norm-1
    
    % update R,S,T
    R = R+mu*(W-H*X);
    S = S+mu*(V-X*OP);
    PX = zeros(size(Z));
    for i=1:F
        PX(:,i)=P(:,:,i)*X(:,i);
    end
    T = T + mu*(PX-Z);

    % update mu
    mu = min(mu*rho,mumax);

    %stop criterion
    object_k = object_kp;
    object_kp = norm(H*X,1)+lambda*norm(X*OP,1);
    
    obj(:,iter) = object_kp;
    stop_c1 = object_k-object_kp;
    stop_c2 = norm(object_k-object_kp)/norm(object_kp);
    if stop_c2<tol && stop_c1>0 && iter>10
        conv = true;
        fprintf('TV3D solver converged at its %4d iterations.\n',iter); 
    elseif iter>maxiter
        conv = true;
        fprintf('TV3D solver reached the maximum %4d iterations.\n',maxiter);
    end        
end

end