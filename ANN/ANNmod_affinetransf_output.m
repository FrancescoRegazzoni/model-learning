function [w,theta] = ANNmod_affinetransf_output(numn,w,theta,A,B,idx)
% Returns the modified ANN, where the idx-th output of the network is an
% affine transformation of the corresponding output of the oringinal ANN:
%                   y_new = A * y_old + B
%
% If idx is not specified, the transformation is applied to all the outputs
%
% NB: works only if BetaOutput=1

if nargin == 5
    idxLastLayer = length(w)-numn(end-1)*numn(end)+1:length(w);
    w(idxLastLayer) = w(idxLastLayer)*A;
    idxLastTheta = length(theta)-numn(end)+1:length(theta);
    theta(idxLastTheta) = theta(idxLastTheta)*A - B;
elseif nargin == 6
    idxLastLayer = length(w)-numn(end-1)*numn(end)+1:length(w);
    W = reshape(w(idxLastLayer),numn(end),numn(end-1));
    W(idx,:) = W(idx,:)*A;
    w(idxLastLayer) = reshape(W,[],1);
    idxLastTheta = length(theta)-numn(end)+idx;
    theta(idxLastTheta) = theta(idxLastTheta)*A - B;
end