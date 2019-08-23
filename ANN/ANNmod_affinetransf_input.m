function [w,theta] = ANNmod_affinetransf_input(numn,w,theta,A,B,idx)
% Returns the modified ANN, where the idx-th input of the network is an
% affine transformation of the corresponding input of the oringinal ANN:
%                   x_new = A * x_old + B
%
% If idx is not specified, the transformation is applied to all the inputs

if nargin == 5
    idxValues = 1:numn(1);
elseif nargin == 6
    idxValues = idx;
end

idxFirstLayer = 1:numn(1)*numn(2);
idxFirstTheta = 1:numn(2);
W = reshape(w(idxFirstLayer),numn(2),numn(1));
for idx = idxValues
    theta(idxFirstTheta) = theta(idxFirstTheta) + W(:,idx)*B/A;
    W(:,idx) = W(:,idx)/A;
end
w(idxFirstLayer) = reshape(W,[],1);