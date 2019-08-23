function [doutdinput,doutdw,doutdtheta,ddoutdinputdw,ddoutdinputdtheta] = ANNSensitivityBis(numn,w,df,ddf,alpha,beta,BetaOutput)

if nargin == 5
    BetaOutput = 0;
end

numw = numn(1:end-1).*numn(2:end);
numwcum = [0;cumsum(numw)];
numncum = [0;cumsum(numn)];
numl = length(numn);

% computation of derivatives wrt alpha, beta
doutdalpha = zeros(numn(end),numncum(end));
doutdbeta = zeros(numn(end),numncum(end));
if BetaOutput
    doutdbeta(:,numncum(end-1)+1:numncum(end)) = eye(numn(end));
else
    doutdalpha(:,numncum(end-1)+1:numncum(end)) = eye(numn(end));
end

for i = numl-1:-1:1    
    currIdx = numncum(i+1)+1:numncum(i+2);
    precIdx = numncum(i)+1:numncum(i+1);
    if ~(BetaOutput && i==(numl-1))
        doutdbeta(:,currIdx) = doutdalpha(:,currIdx) .* df(beta(currIdx))';
        doutdbeta(:,currIdx) = doutdalpha(:,currIdx) .* reshape(df(beta(currIdx)),[1 numn(i+1)]);
    end    
    doutdalpha(:,precIdx) = doutdbeta(:,currIdx)*reshape(w(numwcum(i)+1:numwcum(i+1)),numn(i+1),numn(i));
end

doutdinput = doutdalpha(:,1:numn(1));

% computation of derivatives wrt w, theta
doutdtheta = - doutdbeta(:,numn(1)+1:end);

d sigma_r / d w_ij = alpha_i * f'(beta_j) * d sigma_r / d alpha_j
doutdw = zeros(numn(end),numwcum(end));
idwWij = 1;
for iLayer = 1:numl-1
    for iMittente = 1:numn(iLayer)
        i = numncum(iLayer)+iMittente;
        for iRecettore = 1:numn(iLayer+1)
            j = numncum(iLayer+1)+iRecettore;
            doutdw(:,idwWij) = alpha(i) * doutdbeta(:,j);
            idwWij = idwWij+1;
        end
    end
end
% computation of mixed derivatives
ddoutdalphadw = zeros(numn(end),numncum(end),numwcum(end));
ddoutdbetadw = zeros(numn(end),numncum(end),numwcum(end));
ddoutdalphadtheta = zeros(numn(end),numncum(end),numncum(end)-numn(1));
ddoutdbetadtheta = zeros(numn(end),numncum(end),numncum(end)-numn(1));

for i = numl-1:-1:1
    currIdx = numncum(i+1)+1:numncum(i+2);
    precIdx = numncum(i)+1:numncum(i+1);
    if ~(BetaOutput && i==(numl-1))       
        dfbetacurr = reshape(df(beta(currIdx)),[1 numn(i+1)]);
        
        ddoutdbetadtheta(:,currIdx,:) = ddoutdalphadtheta(:,currIdx,:) .* dfbetacurr;
        ddoutdbetadtheta(:,currIdx,currIdx-numn(1)) = ddoutdbetadtheta(:,currIdx,currIdx-numn(1)) - ...
            doutdalpha(:,currIdx) .* reshape(ddf(beta(currIdx)),[1 numn(i+1)]) .* permute(eye(numn(i+1)),[3,1,2]);    
                        
        ddoutdbetadw(:,currIdx,:) = ddoutdalphadw(:,currIdx,:) .* dfbetacurr;
    end    
    for j = 1:numncum(end)-numn(1)
        ddoutdalphadtheta(:,precIdx,j) = ddoutdbetadtheta(:,currIdx,j)*reshape(w(numwcum(i)+1:numwcum(i+1)),numn(i+1),numn(i));
    end
end

ddoutdinputdtheta = ddoutdalphadtheta(:,1:numn(1),:);
ddoutdinputdw = ddoutdinputdtheta .* reshape(alpha(1:numn(1)),[1 numn(1)]);

ddoutdinputdw = zeros(numn(end),numn(1),numwcum(end));

idwWij = 1;
for iLayer = 1:numl-1
    for iMittente = 1:numn(iLayer)
        i = numncum(iLayer)+iMittente;
        for iRecettore = 1:numn(iLayer+1)
            j = numncum(iLayer+1)+iRecettore;
            ddoutdinputdw(:,:,idwWij) = - alpha(i) * ddoutdinputdtheta(:,:,j-numn(1));            
            if iLayer == 1                
                ddoutdinputdw(:,i,idwWij) = ddoutdinputdw(:,i,idwWij) + doutdbeta(:,j);
            end
            idwWij = idwWij+1;
        end
    end
end

idwWij = 1;
for iMittente = 1:numn(1)
    i = iMittente;
    for iRecettore = 1:numn(2)
        j = numn(1)+iRecettore;
        ddoutdinputdw(:,:,idwWij) = ddoutdinputdtheta(:,:,j) * alpha(i);
        ddoutdinputdw(:,i,idwWij) = ddoutdinputdw(:,i,idwWij) + doutdbeta(:,j);
        idwWij = idwWij+1;
    end
end