function [alpha,phiZero,dphiZero] = LinearSearchScalar(evaluateAlpha,opt,phiZero,dphiZero)

if nargin < 4
    [phiZero,dphiZero] = evaluateAlpha(0);
end

if ~isfield(opt,'dphiZeroThreshold')
    opt.dphiZeroThreshold = -1;
end

fprintf('   LS (alg. %d) --- J_0 = %e --- dphi_0 = %1.5e\n',opt.algorithm,phiZero,dphiZero)

if abs(dphiZero/phiZero) < opt.dphiZeroThreshold
    alpha = 0;
    return
end

if opt.algorithm == 2
    LazyGradientComputation = 1;
    c1=1e-4;
    c2=.9;
    alphamin = 1e-4;
    alphamax = 1e3;
    if isfield(opt,'alpha0')
        alpha=opt.alpha0;
    else
        alpha=1;
    end
    alpha = min(alphamax,max(alphamin,alpha));
    maxIter = 1e1;
    iter = 0;

    alphaOld=0; 
    phiOld=phiZero;
    StopLineSearch=0;
    zoom=0;
    while (~StopLineSearch && iter<maxIter)
        iter = iter+1;
        if LazyGradientComputation
            [phi] = evaluateAlpha(alpha);
        else
            [phi,dphiCurr] = evaluateAlpha(alpha);
        end
        
        fprintf('      LS step %d --- alpha = %1.5e --- J = %1.5e\n',iter,alpha,phi)
        
        if (phi>phiZero+c1*alpha*dphiZero || (phi>= phiOld && alphaOld>0))
            fprintf('      passing to zoom... A(%d,%d)\n',phi>phiZero+c1*alpha*dphiZero,(phi>= phiOld && alphaOld>0))
            alphaLow=alphaOld;
            alphaHigh=alpha;
            phiAlphaLow=phiOld;
            zoom=1;
            StopLineSearch=1;
        else
            if LazyGradientComputation
                [phi,dphiCurr] = evaluateAlpha(alpha);
            end
            if (abs(dphiCurr) <= -c2*dphiZero)
                fprintf('      stopping...\n')
                StopLineSearch=1;
            elseif (dphiCurr>=0)
                fprintf('      passing to zoom... B\n')
                alphaLow=alpha;
                alphaHigh=alphaOld;
                phiAlphaLow=phi;
                zoom=1;
                StopLineSearch=1;
            else
                alphaOld=alpha;
                phiOld=phi;
                alpha=alphaOld*2;
            end
        end
    end
    while (zoom && iter<maxIter)
        iter = iter+1;
        alpha=(alphaLow+alphaHigh)/2.;
        if LazyGradientComputation
            [phi] = evaluateAlpha(alpha);
        else
            [phi,dphiCurr] = evaluateAlpha(alpha);
        end
        
        fprintf('      LS zoom %d --- alpha = %1.5e --- J = %e\n',iter,alpha,phi)
        
        if (phi>phiZero+c1*alpha*dphiZero || phi>= phiAlphaLow)
            alphaHigh=alpha;
        else
            if LazyGradientComputation
                [phi,dphiCurr] = evaluateAlpha(alpha);
            end
            if (abs(dphiCurr) <= -c2*dphiZero)
                fprintf('      stopping...\n')
                zoom=0;
            elseif (dphiCurr*(alphaHigh-alphaLow) >=0)
                alphaHigh=alphaLow;
            end
            alphaLow=alpha;
            phiAlphaLow=phi;
        end
    end
elseif opt.algorithm == 1
    delta = 0.5;
    gamma = 1e-4;
    alpha = 1;
    maxIter = 1e2;
    iter = 0;

    goon = 1;
    while (goon && iter<maxIter)
        iter = iter+1;
%         [phi,dphiCurr] = evaluateAlpha(alpha);
%         if (phi<=phiZero+gamma*alpha*dphiCurr)            
        [phi] = evaluateAlpha(alpha);
        
        fprintf('      LS step %d --- alpha = %1.5e --- J = %1.5e\n',iter,alpha,phi)
        
        if (phi<=phiZero+gamma*alpha*dphiZero)
            goon = 0;
        else
            alpha = alpha*delta;
        end    
    end
else
    error('wrong algorithm index');
end