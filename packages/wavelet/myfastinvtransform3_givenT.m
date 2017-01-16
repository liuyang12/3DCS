function theta = myfastinvtransform3_givenT(w,T_row,T_col,T_t)

% Xin Yuan, 3D dct or wavelet transformation
% initial date: 04/30/2013
[row col T] = size(w);
% if strcmp(spbasis.space,'dct')
%   T_row = dct(eye(row));
%   T_col = dct(eye(col));  
% else
%  level = 3;
%  qmf   = MakeONFilter('Daubechies',8); 
%  sig_level_row = log2(row); 
%  sig_level_col = log2(col); 
%  T_row = get_waveletMatrix(qmf,sig_level_row,level,level);
%  T_col = get_waveletMatrix(qmf,sig_level_col,level,level);
% end
% 
% if strcmp(spbasis.time,'dct')
%     T_t = dct(eye(T));
% else
%     T = 2^(ceil(log2(T)));
% %    level = 2;
% %    qmf   = MakeONFilter('Daubechies',8);
% %    sig_level_t = log2(row);
% %    T_t = get_waveletMatrix(qmf,sig_level_t,level,level);
%  level = 1;
%    qmf   = MakeONFilter('Haar',4);  % Here we use the 'Haar' wavelet and the second parameter now no meaning
%    sig_level_t = ceil(log2(T));
%    T_t = get_waveletMatrix(qmf,sig_level_t,level,level);
% end

theta = shiftdim(reshape(T_t'*reshape(shiftdim(w,2),T,row*col),T,row,col),1);
theta = reshape(T_row'*reshape(theta,row,col*T),row,col,T);
theta = shiftdim(reshape(T_col'*reshape(shiftdim(theta,1),col,row*T),col,T,row),2);


end