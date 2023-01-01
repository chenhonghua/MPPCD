function [s,U,U2] = svdCov(nnIdx, idx, Data, Seed, bi_normals)

%svd
nnPt = Data(:,nnIdx);
C = matrixCompute(nnPt,Seed(:,idx));
[U,S,~] = svd(C);
s = diag(S)/sum(diag(S));

%compute the rotate matrix between svd_n and bilateral_normal
bi_normals(idx,:) = bi_normals(idx,:)/norm(bi_normals(idx,:));
U(:,3) = sign(dot(U(:,3),-Seed(:,idx)))*U(:,3)*(-1);
rotM = vrrotvec2mat(vrrotvec(U(:,3),bi_normals(idx,:)'));
U2 = rotM*U;



% n1 = sign(dot(U(:,3),-Seed(:,idx)))*U(:,3)*(-1);
% n2 = bi_normals(idx,:)';
% %rotation angle
% rotThet = acos(dot(n1, n2));
% rotAxis = cross(n1,n2);
% rotAxis = rotAxis/norm(rotAxis);
% rotM(1,1) = cos(rotThet) + rotAxis(1)*rotAxis(1)*(1-cos(rotThet));
% rotM(1,2) = rotAxis(1)*rotAxis(2)*(1-cos(rotThet)) - rotAxis(3)*sin(rotThet);
% rotM(1,3) = rotAxis(2)*sin(rotThet) + rotAxis(1)*rotAxis(3)*(1-cos(rotThet));
% 
% rotM(2,1) = rotAxis(3)*sin(rotThet) + rotAxis(1)*rotAxis(2)*(1-cos(rotThet));
% rotM(2,2) = cos(rotThet) + rotAxis(2)*rotAxis(2)*(1-cos(rotThet));
% rotM(2,3) = rotAxis(2)*rotAxis(3)*(1-cos(rotThet)) - rotAxis(1)*sin(rotThet);
% 
% rotM(3,1) = rotAxis(1)*rotAxis(3)*(1-cos(rotThet)) - rotAxis(2)*sin(rotThet) ;
% rotM(3,2) = rotAxis(1)*sin(rotThet) + rotAxis(2)*rotAxis(3)*(1-cos(rotThet));
% rotM(3,3) = cos(rotThet) + rotAxis(3)*rotAxis(3)*(1-cos(rotThet));






