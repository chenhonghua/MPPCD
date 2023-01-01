function [ nnPT_tanslate_rotate, grids, gridsID] = localPara( nnIdx, idx, frame, Data, Seed, gridSize)
%LOCALPARA apply the local frame and get the height map
%   此处显示详细说明
nnPt = Data(:,nnIdx);

%translate to the center(center is the seed point)
nnPT_tanslate = bsxfun(@minus, nnPt, Seed(:,idx));

%rotate
nnPT_tanslate_rotate = frame'*nnPT_tanslate; %note the frame' or frame
nnPT_tanslate_rotate = nnPT_tanslate_rotate';

%init grid
grids = zeros(gridSize*gridSize,1);
gridsID = zeros(gridSize*gridSize,6); % 6 is the number of k nearest neighbour

%for each grid center, we compute the average height of its neighbors,we
%can also use the outlier-robust l1-norm instead of l2-norm
%note: grid is created by the bounding rectangular 
minX = min(nnPT_tanslate_rotate(:,1)); maxX = max(nnPT_tanslate_rotate(:,1));
minY = min(nnPT_tanslate_rotate(:,2)); maxY = max(nnPT_tanslate_rotate(:,2));
sideLenX = maxX - minX;
sideLenY = maxY - minY;
for i = 0 : gridSize*gridSize-1
    ii = floor(i/gridSize) + 1;
    jj = mod(i, gridSize) +1;
    centerX = minX+(ii-1)*sideLenX/(gridSize-1);
    centerY = maxY-(jj-1)*sideLenY/(gridSize-1);
    center = [centerX centerY 0];
    %neighIdx_local = rangesearch(nnPT_tanslate_rotate, center, 1.5);
    [n,d] = knnsearch(nnPT_tanslate_rotate, center, 'k', 6); %average nearest k point heights as the grid value
    grids(i+1,1) = sum((nnPT_tanslate_rotate(n,3)))/6;%sum(abs(nnPT_tanslate_rotate(n,3)))/6;%
    gridsID(i+1,:) = nnIdx(n);
end

% frame2 = frame;
% frame2(:,1) = frame2(:,1)*(-1);
% nnPT_tanslate_rotate = frame2'*nnPT_tanslate; %note the frame' or frame
% nnPT_tanslate_rotate = nnPT_tanslate_rotate';
% minX = min(nnPT_tanslate_rotate(:,1)); maxX = max(nnPT_tanslate_rotate(:,1));
% minY = min(nnPT_tanslate_rotate(:,2)); maxY = max(nnPT_tanslate_rotate(:,2));
% sideLenX = maxX - minX;
% sideLenY = maxY - minY;
% for i = 0 : gridSize*gridSize-1
%     ii = floor(i/gridSize) + 1;
%     jj = mod(i, gridSize) +1;
%     centerX = minX+(ii-1)*sideLenX/(gridSize-1);
%     centerY = maxY-(jj-1)*sideLenY/(gridSize-1);
%     center = [centerX centerY 0];
%     %neighIdx_local = rangesearch(nnPT_tanslate_rotate, center, 1.5);
%     [n,d] = knnsearch(nnPT_tanslate_rotate, center, 'k', 6); %average nearest k point heights as the grid value
%     grids(i+1+gridSize*gridSize*1,1) = sum(abs(nnPT_tanslate_rotate(n,3)))/6;
%     gridsID(i+1,:) = nnIdx(n);
% end
% 
% frame3 = frame;
% frame3(:,2) = frame3(:,2)*(-1);
% nnPT_tanslate_rotate = frame3'*nnPT_tanslate; %note the frame' or frame
% nnPT_tanslate_rotate = nnPT_tanslate_rotate';
% minX = min(nnPT_tanslate_rotate(:,1)); maxX = max(nnPT_tanslate_rotate(:,1));
% minY = min(nnPT_tanslate_rotate(:,2)); maxY = max(nnPT_tanslate_rotate(:,2));
% sideLenX = maxX - minX;
% sideLenY = maxY - minY;
% for i = 0 : gridSize*gridSize-1
%     ii = floor(i/gridSize) + 1;
%     jj = mod(i, gridSize) +1;
%     centerX = minX+(ii-1)*sideLenX/(gridSize-1);
%     centerY = maxY-(jj-1)*sideLenY/(gridSize-1);
%     center = [centerX centerY 0];
%     %neighIdx_local = rangesearch(nnPT_tanslate_rotate, center, 1.5);
%     [n,d] = knnsearch(nnPT_tanslate_rotate, center, 'k', 6); %average nearest k point heights as the grid value
%     grids(i+1+gridSize*gridSize*2,1) = sum(abs(nnPT_tanslate_rotate(n,3)))/6;
%     gridsID(i+1,:) = nnIdx(n);
% end
% 
% frame4 = frame;
% frame4(:,1:2) = frame4(:,1:2)*(-1);
% nnPT_tanslate_rotate = frame4'*nnPT_tanslate; %note the frame' or frame
% nnPT_tanslate_rotate = nnPT_tanslate_rotate';
% minX = min(nnPT_tanslate_rotate(:,1)); maxX = max(nnPT_tanslate_rotate(:,1));
% minY = min(nnPT_tanslate_rotate(:,2)); maxY = max(nnPT_tanslate_rotate(:,2));
% sideLenX = maxX - minX;
% sideLenY = maxY - minY;
% for i = 0 : gridSize*gridSize-1
%     ii = floor(i/gridSize) + 1;
%     jj = mod(i, gridSize) +1;
%     centerX = minX+(ii-1)*sideLenX/(gridSize-1);
%     centerY = maxY-(jj-1)*sideLenY/(gridSize-1);
%     center = [centerX centerY 0];
%     %neighIdx_local = rangesearch(nnPT_tanslate_rotate, center, 1.5);
%     [n,d] = knnsearch(nnPT_tanslate_rotate, center, 'k', 6); %average nearest k point heights as the grid value
%     grids(i+1+gridSize*gridSize*3,1) = sum(abs(nnPT_tanslate_rotate(n,3)))/6;
%     gridsID(i+1,:) = nnIdx(n);
% end

%old
% for i = 1 : gridSize
%     for j = 1 : gridSize
%         centerX = minX+(i-1)*sideLenX/(gridSize-1);
%         centerY = maxY-(j-1)*sideLenY/(gridSize-1);
%         center = [centerX centerY 0];
%         %neighIdx_local = rangesearch(nnPT_tanslate_rotate, center, 1.5);
%         [n,d]=knnsearch(nnPT_tanslate_rotate, center, 'k', 1);
%         grids(i,j) = d;
%         gridsID(i,j) = nnIdx(n);
%     end
% end
end

