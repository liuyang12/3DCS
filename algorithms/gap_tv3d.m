function  [ x_3d ] = gap_tv3d( y,A_fun,At_fun,params )
%GAP_TV3D Total variation (TV)-based generalized alternating projection
%method in three dimensions (3D).
%   x_3d=GAP_TV3D(y,A_fun,At_fun,params) returns the reconstructed three 
%   dimensional volume x_3d, which follows the forward model of three 
%   dimensional compressive sensing (3DCS) scheme i.e. y=A_fun(x_3d), where
%   A_fun is the forward model of three dimensional compressive sensing
%   FM_3DCS, y is the M-by-F measurement matrix, x_3d is the (m-by-n)-by-F 
%   three dimensional signal, and params are the parameters used in the
%   algorithm.
%   [model]
%     min  |TV_3(x_3d)|_1
%     s.t. y=A_fun(x_3d)
%   See also TVDENOISING, CS3D.

% [0] parameters configuration 
lambda   = params.lambda;
tvweight = params.tvweight;
tviter   = params.tviter;
maxiter  = params.maxiter;

rows = params.rows;
cols = params.cols;
nframe = params.nframe;

Phi_sum = A_fun(At_fun(ones(size(y)))); % [cs] replaced by function handle

% [1] initialization
theta_3d = At_fun(y);
yt       = zeros(size(y));
for iter = 1:maxiter
    yb = A_fun(reshape(theta_3d,[rows*cols,nframe]));
    if(params.acc)
        yt = yt + (y-yb);
        theta_3d = theta_3d + lambda.*At_fun((yt-yb)./Phi_sum);
    else
        theta_3d = theta_3d + lambda.*At_fun((y-yb)./Phi_sum);
    end
    % TV denoising-based update of theta
    theta_3d = tvdenoising(theta_3d,tvweight,tviter);  
end

x_3d = theta_3d;

end

