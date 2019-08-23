function [varargout] = ANNeval(ANN,input,outputs,algorithm,alpha,beta)

% tic
numn        = ANN.numn       ;
w           = ANN.w          ;
theta       = ANN.theta      ;
f           = ANN.f          ;
df          = ANN.df         ;
ddf         = ANN.ddf        ;
BetaOutput  = ANN.BetaOutput ;

computeAlphaBeta = (nargin == 4);
if strcmp(algorithm,'F') || strcmp(algorithm,'FP')
   ForwP = 1; 
elseif strcmp(algorithm,'B') || strcmp(algorithm,'BP')
   ForwP = 0; 
else
    error('Unrecognized algorithm')
end

Outa  = any(strcmp(outputs,'a'));
Outb  = any(strcmp(outputs,'b'));
Outi  = any(strcmp(outputs,'i'));
Outw  = any(strcmp(outputs,'w'));
Outt  = any(strcmp(outputs,'t'));
Outiw = any(strcmp(outputs,'iw'));
Outit = any(strcmp(outputs,'it'));
if Outiw
    Outi = 1;
    Outw = 1;
end
if Outit
    Outi = 1;
    Outt = 1;
end
Computedf = Outi || Outw || Outt;
Computeddf = Outiw || Outit;
ComputeFullAlpha = Outa || Outw || ~ForwP;  % il calcolo di dinputdw richiede tutta la matrice delle alpha
ComputeFullBeta = Outa || ~ForwP;

if (Outiw || Outit) && ~ForwP
    error('Mixed derivatives are not implemented with Backward Propagation')
end

numw = numn(1:end-1).*numn(2:end);
numwcum = [0;cumsum(numw)];
numncum = [0;cumsum(numn)];
numl = length(numn);

nN = numncum(end);
nI = numn(1);
nW = numwcum(end);
nT = nN-nI;

if computeAlphaBeta 
    if ComputeFullAlpha alpha = zeros(numncum(end),1); end
    if ComputeFullBeta beta = zeros(numncum(end),1); end
end

if computeAlphaBeta
    alphaCurr = input;
    if ComputeFullAlpha
        alpha(1:numn(1)) = input;
    end
end

if ForwP
    if Outi
        dalphadinput = zeros(nN,nI);
        dbetadinput = zeros(nN,nI);
        dalphadinput(1:nI,1:nI) = eye(nI);
    end
    if Outw
        dalphadw = zeros(nN,nW);
        dbetadw = zeros(nN,nW);
    end
    if Outt
        dalphadtheta = zeros(nN,nT);
        dbetadtheta = zeros(nN,nT);
        dbetadtheta(nI+1:end,:) = -eye(nT);
    end
    if Outiw
        ddalphadinputdw = zeros(nN,nI,nW);
        ddbetadinputdw = zeros(nN,nI,nW);
    end
    if Outit
        ddalphadinputdtheta = zeros(nN,nI,nT);
        ddbetadinputdtheta = zeros(nN,nI,nT);
    end
    idwWij = 1;
end
% T1=toc
% tic
if computeAlphaBeta || ForwP
    for iLayer = 2:numl 
        currIdx = numncum(iLayer)+1:numncum(iLayer+1);
        precIdx = numncum(iLayer-1)+1:numncum(iLayer);
        Wcurr = reshape(w(numwcum(iLayer-1)+1:numwcum(iLayer)),numn(iLayer),numn(iLayer-1));

        if computeAlphaBeta 
            betaCurr = Wcurr*alphaCurr - theta(currIdx-numn(1));
            alphaCurr = f(betaCurr);
            if ComputeFullAlpha
                alpha(currIdx) = alphaCurr;
            end
            if ComputeFullBeta        
                beta(currIdx) = betaCurr;
            end
        else
            alphaCurr = alpha(currIdx);
            betaCurr = beta(currIdx);
        end

        if ForwP
            if Computedf
                dfbetacurr = df(betaCurr);
            end
            if Computeddf
                ddfbetacurr = ddf(betaCurr);
            end

            if Outi
                dbetadinput(currIdx,:) = Wcurr*dalphadinput(precIdx,:);
                dalphadinput(currIdx,:) = dbetadinput(currIdx,:).*dfbetacurr;
            end
            if Outw
                dbetadw(currIdx,:) = Wcurr*dalphadw(precIdx,:);
            end
            if Outiw
                ddbetadinputdw(currIdx,:,:) = reshape(Wcurr*reshape(ddalphadinputdw(precIdx,:,:),numn(iLayer-1),[]),numn(iLayer),nI,nW);
            end

            if Outw
                for iMittente = 1:numn(iLayer-1)
                    i = numncum(iLayer-1)+iMittente;
                    for iRecettore = 1:numn(iLayer)
                        j = numncum(iLayer)+iRecettore;
                        dbetadw(j,idwWij) = dbetadw(j,idwWij) + alpha(i);
                        if Outiw
                            ddbetadinputdw(j,:,idwWij) = ddbetadinputdw(j,:,idwWij) + dalphadinput(i,:);
                        end
                        idwWij = idwWij+1;
                    end
                end
            end

            if Outw
                dalphadw(currIdx,:) = dbetadw(currIdx,:).*dfbetacurr;    
            end
            if Outiw
                ddalphadinputdw(currIdx,:,:) = ddfbetacurr.*dbetadinput(currIdx,:).*permute(dbetadw(currIdx,:),[1 3 2]) ...
                    + dfbetacurr.*ddbetadinputdw(currIdx,:,:);
            end    
            if Outt
                dbetadtheta(currIdx,:) = dbetadtheta(currIdx,:) + Wcurr*dalphadtheta(precIdx,:);
                dalphadtheta(currIdx,:) = dbetadtheta(currIdx,:).*dfbetacurr;
            end
            if Outit
                ddbetadinputdtheta(currIdx,:,:) = reshape(Wcurr*reshape(ddalphadinputdtheta(precIdx,:,:),numn(iLayer-1),[]),numn(iLayer),nI,nT);
                ddalphadinputdtheta(currIdx,:,:) = ddfbetacurr.*dbetadinput(currIdx,:).*permute(dbetadtheta(currIdx,:),[1 3 2]) ...
                    + dfbetacurr.*ddbetadinputdtheta(currIdx,:,:);
            end
        end
    end
end

if ~ForwP
    if Outi || Outw || Outt        
        doutdalpha = zeros(numn(end),numncum(end));
        if ~BetaOutput
            doutdalpha(:,numncum(end-1)+1:numncum(end)) = eye(numn(end));
        end

        for i = numl-1:-1:1
            if BetaOutput && i==(numl-1)
                doutdalpha(:,numncum(i)+1:numncum(i+1)) = reshape(w(numwcum(i)+1:numwcum(i+1)),numn(i+1),numn(i));
            else
                doutdalpha(:,numncum(i)+1:numncum(i+1)) = doutdalpha(:,numncum(i+1)+1:numncum(i+2)) * ...
                    (reshape(w(numwcum(i)+1:numwcum(i+1)),numn(i+1),numn(i)) .* ...
                    repmat(df(beta(numncum(i+1)+1:numncum(i+2))),1,numn(i)) );
            end
        end
    end
    
    if Outt
        doutdtheta = - repmat(df(beta(numn(1)+1:end))',numn(end),1) .* doutdalpha(:,numn(1)+1:end);
        doutdtheta(:,end-numn(end)+1:end) = - eye(numn(end));
    end
    
    if Outw
        %d sigma_r / d w_ij = alpha_i * f'(beta_j) * d sigma_r / d alpha_j
        idwWij = 1;
        for iLayer = 1:numl-1
            for iMittente = 1:numn(iLayer)
                for iRecettore = 1:numn(iLayer+1)
                    i = numncum(iLayer)+iMittente;
                    j = numncum(iLayer+1)+iRecettore;
                    if BetaOutput && iLayer==(numl-1)
                        doutdw(iRecettore,idwWij) = alpha(i);
                    else
                        doutdw(:,idwWij) = alpha(i) * df(beta(j)) * doutdalpha(:,j);
                    end
                    idwWij = idwWij+1;
                end
            end
        end
    end
end
% T2=toc
% tic
for iOut=1:length(outputs)
    switch outputs{iOut}
        case 'a'
            varargout{iOut} = alpha; 
        case 'b'
            varargout{iOut} = beta;
        case 'o'        
            if ForwP
                if BetaOutput
                    varargout{iOut} = betaCurr;
                else
                    varargout{iOut} = alphaCurr;
                end 
            else
                if BetaOutput
                    varargout{iOut} = beta(numncum(numl-1)+1:numncum(numl));
                else
                    varargout{iOut} = alpha(numncum(numl-1)+1:numncum(numl));
                end 
            end
        case 'i'         
            if ForwP
                if BetaOutput
                    varargout{iOut} = dbetadinput(currIdx,:);
                else
                    varargout{iOut} = dalphadinput(currIdx,:);
                end 
            else
                varargout{iOut} = doutdalpha(:,1:numn(1));
            end
        case 'w'                
            if ForwP
                if BetaOutput
                    varargout{iOut} = dbetadw(currIdx,:);
                else
                    varargout{iOut} = dalphadw(currIdx,:);
                end 
            else
                varargout{iOut} = doutdw;
            end
        case 't'         
            if ForwP
                if BetaOutput
                    varargout{iOut} = dbetadtheta(currIdx,:);
                else
                    varargout{iOut} = dalphadtheta(currIdx,:);
                end 
            else
                varargout{iOut} = doutdtheta;
            end
        case 'iw'
            if BetaOutput
                varargout{iOut} = ddbetadinputdw(currIdx,:,:);
            else
                varargout{iOut} = ddalphadinputdw(currIdx,:,:);
            end  
        case 'it'
            if BetaOutput
                varargout{iOut} = ddbetadinputdtheta(currIdx,:,:);
            else
                varargout{iOut} = ddalphadinputdtheta(currIdx,:,:);
            end 
    end
end
% T3=toc
% 
% T1/T2
% T3/T2