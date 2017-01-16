function weight = my_get_weight3d_tree(row,col,T,block)

% Xin Yuan, get the weight of coefficient for GAP 3D
% initial date: 04/30/2013
% Revised date: 05/04/2013, add the zero-tree structure of the weight for
% wavelet tree


% Input: row col, T the video size or the coeffient size
% Input: block, the blocksize

% output: weight, will have the size [block_num_row block_num_col block_num_T], in order to weight the
% coefficient one by one
% Here we use the cosine weight to weight the coefficient the low frequency
% will have the larger weight


row_block = block.row;
col_block = block.col;
T_block = block.T;

block_num_row = length(row_block);
block_num_T = ceil(T/T_block);
%block_num_col = length(col_block);
weight = zeros(row,col,T);

weight_2D = zeros(row,col);
weight_row = 1e-5+cos((0:(block_num_row-1))/block_num_row*pi/2).^2;
weight_T = 1e-5+cos((0:(block_num_T-1))/block_num_T*pi/2).^2;
%weight_col = 1e-5+cos((0:(block_num_col-1))/block_num_col*pi/2).^2;

weight_2D((1:row_block(1)), (1:col_block(1))) = kron(weight_row(1),ones(row_block(1),col_block(1)));

rowsum = row_block(1);
colsum = col_block(1);
for nR = 2:block_num_row
    
    temp_weight =  kron(weight_row(nR), ones(row_block(nR),col_block(nR)));
    weight_2D(rowsum+(1:row_block(nR)), colsum+(1:col_block(nR))) = temp_weight;
    weight_2D(rowsum+(1:row_block(nR)), 1:colsum) = temp_weight;
    weight_2D(1:rowsum, colsum+(1:col_block(nR))) = temp_weight;
    
    rowsum = rowsum + row_block(nR);
    colsum = colsum + col_block(nR);
end

for nt = 1:block_num_T
    if(nt*T_block<=T)
        weight(:,:,(nt-1)*T_block+(1:T_block))  = repmat(weight_2D*weight_T(nt), [1 1 T_block]);
    else
        weight(:,:,((nt-1)*T_block+1):T) = repmat(weight_2D*weight_T(nt), [1, 1, T-((nt-1)*T_block)]);
    end
end


end