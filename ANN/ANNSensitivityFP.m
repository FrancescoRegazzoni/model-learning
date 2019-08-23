function [doutdinput,doutdw,doutdtheta,ddoutdinputdw,ddoutdinputdtheta] = ANNSensitivityFP(numn,w,df,ddf,alpha,beta,BetaOutput)

if nargin == 5
    BetaOutput = 0;
end

numw = numn(1:end-1).*numn(2:end);
numwcum = [0;cumsum(numw)];
numncum = [0;cumsum(numn)];
numl = length(numn);

nN = numncum(end);
nO = numn(end);
nI = numn(1);
nW = numwcum(end);
nT = nN-nI;

dalphadinput = zeros(nN,nI);
dbetadinput = zeros(nN,nI);
dalphadw = zeros(nN,nW);
dbetadw = zeros(nN,nW);
dalphadtheta = zeros(nN,nT);
dbetadtheta = zeros(nN,nT);
dbetadtheta(nI+1:end,:) = -eye(nT);
ddalphadinputdw = zeros(nN,nI,nW);
ddbetadinputdw = zeros(nN,nI,nW);
ddalphadinputdtheta = zeros(nN,nI,nT);
ddbetadinputdtheta = zeros(nN,nI,nT);

dalphadinput(1:nI,1:nI) = eye(nI);

idwWij = 1;
for iLayer = 2:numl   
    currIdx = numncum(iLayer)+1:numncum(iLayer+1);
    precIdx = numncum(iLayer-1)+1:numncum(iLayer);
    Wcurr = reshape(w(numwcum(iLayer-1)+1:numwcum(iLayer)),numn(iLayer),numn(iLayer-1));
    
    dfbetacurr = df(beta(currIdx));
    ddfbetacurr = ddf(beta(currIdx));
    
    dbetadinput(currIdx,:) = Wcurr*dalphadinput(precIdx,:);
    dalphadinput(currIdx,:) = dbetadinput(currIdx,:).*dfbetacurr;
    
    dbetadw(currIdx,:) = Wcurr*dalphadw(precIdx,:);
    ddbetadinputdw(currIdx,:,:) = reshape(Wcurr*reshape(ddalphadinputdw(precIdx,:,:),numn(iLayer-1),[]),numn(iLayer),nI,nW);
    
    for iMittente = 1:numn(iLayer-1)
        i = numncum(iLayer-1)+iMittente;
        for iRecettore = 1:numn(iLayer)
            j = numncum(iLayer)+iRecettore;
            dbetadw(j,idwWij) = dbetadw(j,idwWij) + alpha(i);
            ddbetadinputdw(j,:,idwWij) = ddbetadinputdw(j,:,idwWij) + dalphadinput(i,:);
            idwWij = idwWij+1;
        end
    end
    dalphadw(currIdx,:) = dbetadw(currIdx,:).*dfbetacurr;    
    ddalphadinputdw(currIdx,:,:) = ddfbetacurr.*dbetadinput(currIdx,:).*permute(dbetadw(currIdx,:),[1 3 2]) ...
        + dfbetacurr.*ddbetadinputdw(currIdx,:,:);
    
    dbetadtheta(currIdx,:) = dbetadtheta(currIdx,:) + Wcurr*dalphadtheta(precIdx,:);
    dalphadtheta(currIdx,:) = dbetadtheta(currIdx,:).*dfbetacurr;
        
    ddbetadinputdtheta(currIdx,:,:) = reshape(Wcurr*reshape(ddalphadinputdtheta(precIdx,:,:),numn(iLayer-1),[]),numn(iLayer),nI,nT);
    ddalphadinputdtheta(currIdx,:,:) = ddfbetacurr.*dbetadinput(currIdx,:).*permute(dbetadtheta(currIdx,:),[1 3 2]) ...
        + dfbetacurr.*ddbetadinputdtheta(currIdx,:,:);
    
end

if BetaOutput
    doutdinput = dbetadinput(currIdx,:);
    doutdw = dbetadw(currIdx,:);
    doutdtheta = dbetadtheta(currIdx,:);
    ddoutdinputdw = ddbetadinputdw(currIdx,:,:);
    ddoutdinputdtheta = ddbetadinputdtheta(currIdx,:,:);
else
    doutdinput = dalphadinput(currIdx,:);
    doutdw = dalphadw(currIdx,:);
    doutdtheta = dalphadtheta(currIdx,:);
    ddoutdinputdw = ddalphadinputdw(currIdx,:,:);
    ddoutdinputdtheta = ddalphadinputdtheta(currIdx,:,:);
end
