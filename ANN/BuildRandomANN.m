function [w,theta] = BuildRandomANN(numn,fixed)

if nargin == 1
    fixed = 0;
end

numw = sum(numn(1:end-1).*numn(2:end));
numtheta = sum(numn(2:end));

if fixed
    rnd('default')
end

w = randn(numw,1);
theta = randn(numtheta,1);


% 
% if fixed
%     load('rvec.mat','rvec');
%     w = rvec(1:numw);
%     theta = rvec(numw+1:numw+numtheta,1);
% else
%     w = randn(numw,1);
%     theta = randn(numtheta,1);
%     
%     %w = 10*w;
%     
% %     w(end-numn(end-1)*numn(end):end) = 1e-2*w(end-numn(end-1)*numn(end):end);
% %     theta(end-numn(end):end)         = 1e-2*theta(end-numn(end):end);
% 
%     % w = zeros(numw,1);
%     % theta = zeros(numtheta,1);
% 
%     % w = ones(numw,1);
%     % theta = zeros(numtheta,1);
% 
%     % w = ones(numw,1)*1e-10;
%     % theta = zeros(numtheta,1);
% end