function [ nnPT_rotate_tanslate ] = local2global( nnIdx, idx, frame, pdata, srcSeed_localPara, E_local_grids_height)
%LOCAL2GLOBAL 此处显示有关此函数的摘要
%   此处显示详细说明
nnPt = srcSeed_localPara{idx};

%change height
nnPt(:,3) = E_local_grids_height(nnIdx);

%rotate
nnPT_rotate = inv(frame')*nnPt'; %note the frame' or frame

%translate to the global center
nnPT_rotate_tanslate = bsxfun(@plus, nnPT_rotate, pdata(:,idx));

end

