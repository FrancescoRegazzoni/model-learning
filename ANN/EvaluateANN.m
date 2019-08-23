function [output,alpha,beta] = EvaluateANN(numn,w,theta,input,f,BetaOutput)

if nargin == 5
    BetaOutput = 0;
end

numw = numn(1:end-1).*numn(2:end);
numwcum = [0;cumsum(numw)];
numncum = cumsum(numn);
numl = length(numn);

alpha = zeros(numncum(end),1);
beta = zeros(numncum(end),1);
alpha(1:numn(1)) = input;
alphaPrec = input;
for i = 2:numl
	betaPrec = reshape(w(numwcum(i-1)+1:numwcum(i)),numn(i),numn(i-1))*alphaPrec - theta((numncum(i-1)+1:numncum(i))-numn(1));
    alphaPrec = f(betaPrec);
    beta(numncum(i-1)+1:numncum(i)) = betaPrec;
    alpha(numncum(i-1)+1:numncum(i)) = alphaPrec;
end
if BetaOutput
    output = betaPrec;
else
    output = alphaPrec;
end
