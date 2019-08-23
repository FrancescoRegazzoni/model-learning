function [numnNew,wNew,thetaNew] = ANNmod_extractvariables_input(numn,w,theta,idx)
% Returns a new network where just the variables of index idx are kept in
% the input (the others are assumed to be 0)

idxFirstLayer = 1:numn(1)*numn(2);
W = reshape(w(idxFirstLayer),numn(2),numn(1));
% idxFirstTheta = 1:numn(2);
% th = theta(idxFirstTheta);

% th_new = th(idx);
W_new = W(:,idx);

numnNew = numn;
numnNew(1) = length(idx);
wNew = [reshape(W_new,[],1);w(numn(1)*numn(2)+1:end)];
% thetaNew = [th_new;theta(numn(2)+1:end)];
thetaNew = theta;