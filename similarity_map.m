local_grids_height_matrix000=local_grids_height_matrix;
aaa=local_grids_height_matrix000-repmat(local_grids_height_matrix000(:,82649), 1,110529);
bbb=aaa.*aaa;
ccc = sum(bbb,1);
ddd = sqrt(ccc);
pdata22 = [pdata;ddd]';
pdata22= double(pdata22);
save pdata22.txt -ascii ddd