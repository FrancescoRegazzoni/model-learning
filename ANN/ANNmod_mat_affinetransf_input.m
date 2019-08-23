function [w,theta] = ANNmod_mat_affinetransf_input(numn,w,theta,A,b)
% Returns the modified ANN, where the idx-th input of the network is an
% affine transformation of the corresponding input of the oringinal ANN:
%                   x_new = A x_old + b
%
% x_2 = W * x_old - th = W A^-1(x_new - b) - th = (W A^-1) x_new - (th - W A^-1 b)

idxFirstLayer = 1:numn(1)*numn(2);
W = reshape(w(idxFirstLayer),numn(2),numn(1));
idxFirstTheta = 1:numn(2);
th = theta(idxFirstTheta);

th_new = th + W*(A\b); % i.e. theta+W*inv(A)*b
W_new = W/A; % i.e. W*inv(A)

w(idxFirstLayer) = reshape(W_new,[],1);
theta(idxFirstTheta) = th_new;