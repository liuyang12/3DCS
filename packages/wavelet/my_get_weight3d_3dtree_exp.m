function weight = my_get_weight3d_3dtree_exp(row,col,T,block,weight_base)

% Xin Yuan, get the weight of coefficient for GAP 3D
% initial date: 04/30/2013
% Revised date: 05/04/2013, add the zero-tree structure of the weight for
% wavelet tree
% Revised date: 05/08/2013, Changed to exponential decay



% Input: row col, T the video size or the coeffient size
% Input: block, the blocksize

% output: weight, will have the size [block_num_row block_num_col block_num_T], in order to weight the
% coefficient one by one
% Here we use the cosine weight to weight the coefficient the low frequency
% will have the larger weight


T = 2^(ceil(log2(T)));

if nargin<5
    weight_base.space = 1.5;
    weight_base.T = 1.5;
end

row_block = block.x;
col_block = block.y;
T_block = block.T;

block_num_row = length(row_block);
block_num_T = length(T_block);
block_3d = max([block_num_row, block_num_T ]);
block_3D_small = min([block_num_row, block_num_T ]);

weight = zeros(row,col,T);

base_row = weight_base.space;
base_T = weight_base.T;

base_3d = max([base_row,base_T]);

weight_1d = 1e-5+base_3d.^(-(0:(block_3d-1)));
  
weight((1:row_block(1)), (1:col_block(1)), (1:T_block(1))) = weight_1d(1); 
rowsum = row_block(1);
colsum = col_block(1);
Tsum = T_block(1);

for nR = 2:block_3D_small
    weight(1:rowsum, colsum+(1:col_block(nR)),1:Tsum) = weight_1d(nR);
    weight(rowsum+(1:row_block(nR)), 1: (colsum+col_block(nR)),1:Tsum) = weight_1d(nR);
    weight(1:(rowsum+row_block(nR)), 1: (colsum+col_block(nR)),Tsum + (1:T_block(nR))) = weight_1d(nR);  % The second block
    rowsum = rowsum + row_block(nR);
    colsum = colsum + col_block(nR);
    Tsum = Tsum + T_block(nR);
end

if (block_num_row>block_3D_small)
    for nR = (block_3D_small+1):block_num_row
        weight(rowsum+(1:row_block(nR)),1:colsum) = weight_1d(nR);
        
        weight(1:(rowsum+row_block(nR)),colsum+(1:col_block(nR))) = weight_1d(nR);
        
        rowsum = rowsum + row_block(nR);
        colsum = colsum + col_block(nR);
    end
end

if(block_num_T>block_3D_small)
    for nR = (block_3D_small+1):block_num_T
        
        weight(:,:,Tsum+(1:T_block(nR))) = weight_1d(nR);
        Tsum = Tsum + T_block(nR);
    
    end
end


end