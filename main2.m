clc;clear all;
for iter = 0 : 10
    fileName = strcat('block-ours-normal_iter',num2str(iter));
    fileName2 = strcat(fileName,'.ply');
    
    %% read point cloud
    pc = pcread(fileName2);
    pc_bi = pcread(fileName2);
    bi_normals = pc_bi.Normal;
    pdata = pc.Location';

    %% find R_neighbor
    neighIdx2 = knnsearch(pdata', pdata', 'k', 60);
    neighIdx = mat2cell(neighIdx2, ones(1,pc.Count));
    %neighIdx = rangesearch(pdata', pdata', 0.026);% girl(0.06) model(0.04) boy(0.03) david(0.03) cube2(0.2) 002(0.06) fig13(0.04) qingtongqi2(0.034) buddha(0.03) gargo(4.34) cube(0.13)  welsh-dragon(4) angle(4) block2(3) iron_1noisy(0.04) captain(0.04)

    %% local coordinate frame k
    idxSz = cellfun(@length,neighIdx,'uni',true);
    srcIdx = neighIdx(idxSz>5);
    srcSeed = pdata(:,idxSz>5);
    M = sum(idxSz>5);
    idx = num2cell((1:M)');
    [s,frame_original, frames_filetered] = cellfun(@(x,y)svdCov(x,y,pdata,srcSeed,bi_normals),srcIdx,idx,'uni',false); %n is the local frame

    %% normal covariance
    N_original = cell2mat(frame_original);
    N_original = N_original(:,3);%normal per point
    N_original2 = reshape(N_original, 3, size(N_original,1)/3);
    
    N_filetered = cell2mat(frames_filetered);
    N_filetered = N_filetered(:,3);
    N_filetered2 = reshape(N_filetered, 3, size(N_filetered,1)/3);

    for i = 1 : M
        nidx = neighIdx{i}; %%邻域索引
        nnormal = N_filetered2(:,nidx); %%根�?�索引得到滤波�?�法�?
        mean_nnormal = sum(nnormal,2)/size(nnormal,2); %%平�?�法线（没用到�?
        Cnormal = matrixCompute(nnormal,nnormal(:,1)); %%法线�??方�?
        [U,S,~] = svd(Cnormal); %%svd分解
        %U(:,3) = sign(dot(U(:,3),-nnormal(:,1)))*U(:,3)*(-1);
        rotM = vrrotvec2mat(vrrotvec(U(:,1),bi_normals(i,:)')); %%计算旋转矩阵
        U = rotM*U; %%旋转
        U(:, [1, 3]) = U(:, [3, 1]);
        UU{i} = U;
    end

    N_cov = cell2mat(UU');
    N_cov = N_cov(:,3);
    N_cov2 = reshape(N_cov, 3, size(N_cov,1)/3);
    
    %% local parameterization
    gridSize = 7;
    [srcSeed_localPara,local_grids_height,local_grids_ID] = cellfun(@(x,y,z)localPara(x,y,z,pdata,srcSeed,gridSize),srcIdx,idx,UU','uni',false);
 
    local_grids_height_matrix = cell2mat(local_grids_height);
    local_grids_height_matrix = reshape(local_grids_height_matrix, gridSize*gridSize, M);
   
    %% find similar patches by height map
    patNum = 40;
    [Init_Similar_Index,Init_Similar_Dist] = knnsearch(local_grids_height_matrix', local_grids_height_matrix', 'k', patNum);
    Init_Similar_Index = Init_Similar_Index';
    
   %% LR with graph constraint
    E_local_grids_height_matrix = zeros(size(local_grids_height_matrix));
    W = zeros(size(local_grids_height_matrix));
    for i = 1 : size(Init_Similar_Index,2)
        group_idx = Init_Similar_Index(:,i);
        
        %random permutation
        %group_idx = group_idx(randperm(numel(group_idx)));
         
        group_points = local_grids_height_matrix(:,group_idx);
        M_Temp = repmat(mean( group_points, 2 ),1,patNum);
        group_points = group_points - M_Temp;
        %E_Temp = WSNM( group_points, 3, 0.002, M_Temp, 3); % WNNM Estimation
        E_Temp = LR_ours( group_points, 100, 100); % Our LR Estimation
        E_local_grids_height_matrix(:,group_idx) = E_local_grids_height_matrix(:,group_idx) + E_Temp;
        W(:,group_idx) = W(:,group_idx)+ones(gridSize*gridSize,size(group_idx,1)); %each patch count 
    %     E_all{i} = E_Temp;
    end
    E_local_grids_height_matrix2 = E_local_grids_height_matrix./W; %averge
  
    %% recovery
    %average local z, or we say height for each patch
    count_all = zeros(size(Init_Similar_Index,2),1);
    E_local_grids_height = zeros(size(Init_Similar_Index,2),1);
    for i = 1 : size(Init_Similar_Index,2)
        group_i = Init_Similar_Index(:,i);
        local_patch_height_i = E_local_grids_height_matrix2(:,i);
        local_patch_height_i = repmat(local_patch_height_i, 1, 6);
        local_patch_id_i = local_grids_ID{i};
        local_patch_id_i = local_patch_id_i(:);
        for j = 1 : size(local_patch_id_i,1)*size(local_patch_id_i,2)
             E_local_grids_height(local_patch_id_i(j)) = E_local_grids_height(local_patch_id_i(j)) + local_patch_height_i(j);
             count_all(local_patch_id_i(j)) = count_all(local_patch_id_i(j)) + 1;
        end 
    end
    E_local_grids_height = E_local_grids_height./count_all;

    %use simialer patch to recover the whole model in the global frame
    Data_localPara_global = cellfun(@(x,y,z)local2global(x,y,z,pdata,srcSeed_localPara,E_local_grids_height),srcIdx,idx,UU','uni',false);
    pdata_final = zeros(size(pdata));
    pdata_W = zeros(3,size(pdata,2));
    for i = 1 : size(pdata,2)
        nidx = neighIdx{i};
        pdata_global_i = Data_localPara_global{i};
        %pdata_global_i = pdata(:,nidx); %no p-norm
        %pdata_final(:,nidx) = pdata_global_i; % just value it once, dont average
        
        %simple average
        %pdata_W(:,nidx) = pdata_W(:,nidx)+1;
        %pdata_final(:,nidx) = pdata_final(:,nidx) + pdata_global_i;
        
        %l1-median
        pdata_final(3*pdata_W(:,nidx)+1:3*pdata_W(:,nidx)+3,nidx) =  pdata_global_i;
        pdata_W(:,nidx) = pdata_W(:,nidx)+1;
    end
    %simple average
    %pdata_final = pdata_final./pdata_W;
    
    %/////l1-media way to update the final points
    tic;
    pdata_final2 = zeros(size(pdata));
    for i = 1 : size(pdata_final2,2)
        pdata_final_i = pdata_final(:,i);
        pdata_final_i(find(pdata_final_i==0))=[];
        pdata_final_i = reshape(pdata_final_i, 3, size(pdata_final_i,1)/3);
     
        x0 = sum(pdata_final_i,2)/size(pdata_final_i,2);
        %iteration
        xj = x0;
        for j = 1 : 10
            dist = pdata_final_i-repmat(xj, 1, size(pdata_final_i,2));
            dist2 = sum(dist.^2,1);
            dist2_sqrt = sqrt(sum(dist.^2,1));
            dist2_gauss = exp(-1*dist2./0.05);
            alpha = dist2_gauss./dist2_sqrt;
            alpha2 = repmat(alpha,3,1);
            xj = (pdata_final_i.*alpha2)/sum(alpha);
            xj = sum(xj,2);
        end   
        pdata_final2(:,i) = xj;
    end
    toc;
    %/////l1-media way to update the final points


    %% save
    pdata_final_pc = pointCloud(pdata_final2');
    pdata_final_pc.Normal = bi_normals;
    
    fileName = [];
    fileName = strcat('block-ours-normal_iter',num2str(iter+1));
    %fileName = strcat(fileName,'_outlier');
    pcwrite(pdata_final_pc, fileName);
    %write_to_ply( pdata_final', 'cube_final.ply');
end








