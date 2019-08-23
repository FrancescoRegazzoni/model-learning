function [numnNew,wNew,thetaNew] = ANNmod_extractvariables_output(numn,w,theta,idx)
% Returns a new network where just the variables of index idx are kept in
% the output

idxLastLayer = length(w)-numn(end-1)*numn(end)+1:length(w);
W = reshape(w(idxLastLayer),numn(end),numn(end-1));
idxLastTheta = length(theta)-numn(end)+1:length(theta);
th = theta(idxLastTheta);

th_new = th(idx);
W_new = W(idx,:);

numnNew = numn;
numnNew(end) = length(idx);
wNew = [w(1:length(w)-numn(end-1)*numn(end));reshape(W_new,[],1)];
thetaNew = [theta(1:length(theta)-numn(end));th_new];