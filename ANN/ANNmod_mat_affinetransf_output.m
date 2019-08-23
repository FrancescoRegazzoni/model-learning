function [w,theta] = ANNmod_mat_affinetransf_output(numn,w,theta,A,b)
% Returns the modified ANN, where the idx-th output of the network is an
% affine transformation of the corresponding output of the oringinal ANN:
%                   y_new = A y_old + b
%
% NB: works only if BetaOutput=1
%
% y_new = A y_old + b = A (W x_{n-1} - th) + b = (A W) x_{n-1} - (A th - b)

idxLastLayer = length(w)-numn(end-1)*numn(end)+1:length(w);
W = reshape(w(idxLastLayer),numn(end),numn(end-1));
idxLastTheta = length(theta)-numn(end)+1:length(theta);
th = theta(idxLastTheta);

W_new = A*W;
th_new = A*th - b;

w(idxLastLayer) = reshape(W_new,[],1);
theta(idxLastTheta) = th_new;