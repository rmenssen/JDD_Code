%find the actual minimum. 

[M,I]=min(map(:))
[I_row, I_col] = ind2sub(size(map),I)
Vs(I_row)
kvs(I_col)