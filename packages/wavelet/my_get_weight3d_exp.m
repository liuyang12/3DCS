function weight = my_get_weight3d_exp(row,col,T,block,weight_base)

% Xin Yuan, get the weight of coefficient for GAP 3D
% initial date: 04/30/2013
% revise date: 05/08/2013  Changed to exponential decay

% Input: row col, T the video size or the coeffient size
% Input: block, the blocksize

% output: weight, will have the size [block_num_row block_num_col block_num_T], in order to weight the
% coefficient one by one
% Here we use the exponential decay with the base given by weight_base to weight the coefficient the low frequency
% will have the larger weight


if nargin<5
    weight_base.space = 2;
    weight_base.T = 2;
end


row_block = block.row;
col_block = block.col;
T_block = block.T;
base_row = weight_base.space;
base_col = weight_base.space;
base_T = weight_base.T;

block_num_row = row/row_block;
block_num_col = col/col_block;
block_num_T = ceil(T/T_block);

% weight_row = 1e-5+cos((0:(block_num_row-1))/block_num_row*pi/2).^2;
% weight_col = 1e-5+cos((0:(block_num_col-1))/block_num_col*pi/2).^2;
% weight_T = 1e-5+cos((0:(block_num_T-1))/block_num_T*pi/2).^2;

weight_row = 1e-5+base_row.^(-(0:(block_num_row-1)));
weight_col = 1e-5+base_col.^(-(0:(block_num_col-1))); 
weight_T = 1e-5+base_T.^(-(0:(block_num_T-1)));

weight_row = kron(weight_row, ones(1,row_block));
weight_col = kron(weight_col, ones(1,col_block));
weight_T = kron(weight_T, ones(1,T_block));

weight_row_all = repmat(weight_row', [1, col, T]);
weight_col_all = repmat(weight_col, [row, 1, T]);
weight_T_all = reshape(repmat(weight_T, [row*col, 1]), [row, col, T] );

% weight_row_all = repmat(weight_row', [1, block_num_col, block_num_T]);
% weight_col_all = repmat(weight_col, [block_num_row, 1, block_num_T]);
% weight_T_all = reshape(repmat(weight_T, [block_num_row*block_num_col, 1]), [block_num_row, block_num_col, block_num_T] );

weight = weight_row_all+weight_col_all + weight_T_all;
%weight = weight_block;
%weight = kron(weight_block, ones(row_block,col_block,T_block));


end