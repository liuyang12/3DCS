function w = myfasttransform3_givenT(theta_temp,T_row,T_col,T_t)

% Xin Yuan, 3D dct or wavelet transformation
% initial date: 04/30/2013
% Revised date: 05/09/2013 make the T domain also wavelet

[row, col, T] = size(theta_temp);
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
% %    sig_level_t = log2(T);
% %    T_t = get_waveletMatrix(qmf,sig_level_t,level,level);
%    level = 1;
%    qmf   = MakeONFilter('Haar',4);  % Here we use the 'Haar' wavelet and the second parameter now no meaning
%   % qmf   = MakeONFilter('Daubechies',4); 
%    sig_level_t = ceil(log2(T));
%    T_t = get_waveletMatrix(qmf,sig_level_t,level,level);
% 
% end

w = shiftdim(reshape(T_t*reshape(shiftdim(theta_temp,2),T,row*col),T,row,col),1);
w = reshape(T_row*reshape(w,row,col*T),row,col,T);
w = shiftdim(reshape(T_col*reshape(shiftdim(w,1),col,row*T),col,T,row),2);

end