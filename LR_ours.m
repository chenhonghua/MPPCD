function [ L ] = LR_ours( G, beta, lamda)
%LR_OURS 此处显示有关此函数的摘要
%   此处显示详细说明

%initial U
[U,SigmaY,VV] = svd(G);
% Sigma = zero(size(SigmaY));
% Sigma(1:20,1:20) = SigmaY(1:20,1:20);
% U = U(:,1:15);

%initial W
I = eye(size(G,2));
W = pdist2(G',G');
W = W*(-1);
W = exp(W/0.2);
% W = 1./(W+0.0001);

%solving V
V = U'*G*inv(I+lamda*I+beta*W');

%solving U
[U2,SigmaY2,V2] = svd(V*G'); 
U = V2*U2';

%L
L = U*V;


% |tr(AB)| <= sum(alpha_i*beta_i);
% B=rand(10,10);
% [U,bata,V] = svd(B);
% 
% A = U;
% [U,alpha,V] = svd(B);

end

