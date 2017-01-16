function Obj_rec = run_gap_cassi(y, Phi)
%% use GAP to reconstruct
[ro, co, T] = size(Phi);
rlevel = ceil(log2(ro));
clevel = ceil(log2(co));
rows = 2^rlevel;
cols = 2^clevel;
Phi_new = zeros(rows, cols,T);
Phi_new(1:ro, 1:co,:) = Phi;
clear Phi
Phi = Phi_new; 
clear Phi_new

[Row, Col, T] = size(Phi);
y0 = zeros(rows, cols);
y0(1:ro, 1:co) = y;
clear y
y = y0;
clear y0;
% addpath('./GAP')
spbasis.space = 'wavelet';  % transform for space, 'wavelet' or 'dct'
spbasis.time  = 'dct';  % transform for spectrum, 'wavelet' or 'dct', dct is always used. I f we use wavelet, T need to be the power of 2. we can use no, means no transformation
%spbasis.spectrum  = 'no';  % Here we use no, means no transfromation in spectrum

weighttype.space = 'tree';   % Here we can select:  'tree' or 'block'
weighttype.time = 'block';   % Here we can select:  'tree' or 'block', if we use tree, T need to be the power of 2


weight_base.type = 'cos'; % here we can select 'exp' or 'cos'
if strcmp(weight_base.type,'exp')
    weight_base.space = 1.5;   % This weight is the base of exponential decay. should be larger than 1 [1 2] is always used
    weight_base.time = 1.5;
    weight_base.T = 1.5;
end


ydim = prod(size(y));
% The block size for group
block.row = 2;
block.col = 2;
block.T = T/2;

stopc.iternum = 50;
stopc.err = 10^-5;
acc = 2;
ydim = Row*Col;

m_star = ceil(ydim/(block.row*block.col*block.T));

% function handle of forward model of CASSI
A_fun  = @(f) fm_cassi(Phi,f);
At_fun = @(y) fmt_cassi(Phi,y);

theta0 = GAP_3D_wL21_grayscale(y,A_fun,At_fun,Row,Col,T,block,spbasis,m_star,stopc,acc,weight_base,weighttype); 

% [cs] modified as general cs input
% [Row,Col,T] = size(Phi);
% Phi_cs = zeros(Row*Col,Row*Col,T);
% for t = 1:T
%     Phi_slice = Phi(:,:,t);
%     Phi_cs(:,:,t) = diag(Phi_slice(:));
% end
% A_cs  = @(x) fm_cs(Phi_cs,x);
% At_cs = @(y) fmt_cs(Phi_cs,y);
% y_cs  = y(:);
% 
% theta0 = GAP_3D_wL21_grayscale(y_cs,A_cs,At_cs,Row,Col,T,block,spbasis,m_star,stopc,acc,weight_base,weighttype); 

Obj_rec  = theta0(1:ro, 1:co,:);
end