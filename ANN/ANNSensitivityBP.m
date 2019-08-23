function [doutdalpha,doutdw,doutdtheta] = ANNSensitivityBP(numn,w,df,alpha,beta,BetaOutput)

if nargin == 5
    BetaOutput = 0;
end

 dfbeta = df(beta);
%NB: in caso di funzione di attivazione tanh, siccome alpha = f(beta) si
%possono utilizzare le formule sotto. Tuttavia con l'attuale implementazione
%il tempo di valutazione di df è inifluente sul costo totale.
%dfbeta = 1-alpha.^2;
%dfbeta = 1-alpha.*alpha;

numw = numn(1:end-1).*numn(2:end);
numwcum = [0;cumsum(numw)];
numncum = [0;cumsum(numn)];
numl = length(numn);

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
            repmat(dfbeta(numncum(i+1)+1:numncum(i+2)),1,numn(i)) );
    end
end

doutdtheta = - repmat(dfbeta(numn(1)+1:end)',numn(end),1) .* doutdalpha(:,numn(1)+1:end);
doutdtheta(:,end-numn(end)+1:end) = - eye(numn(end));

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
                doutdw(:,idwWij) = alpha(i) * dfbeta(j) * doutdalpha(:,j);
            end
            idwWij = idwWij+1;
        end
    end
end