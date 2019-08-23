function [xnew,alphaStar,phiZero,dphiZero] = linearsearch_vectorial(x0,direction,FuncGrad,opt,phiZero,dphiZero)

if nargin < 6
    [phiZero,dphiZero] = evaluateAlphaVect(0);
end

if dphiZero > 0
    warning('Increasing direction selected! The direction will be reverted')
    direction = -direction;
    dphiZero = -dphiZero;
end

[alphaStar,phiZero,dphiZero] = LinearSearchScalar(@evaluateAlphaVect,opt,phiZero,dphiZero);
xnew = x0 + alphaStar*direction;    

function [phi,dphi] = evaluateAlphaVect(alphaCurr)
    x = x0 + alphaCurr*direction;    
    if nargout == 1
        [phi] = FuncGrad(x);
    else        
        [phi,DE] = FuncGrad(x);
        dphi = DE'*direction;
    end
end

end