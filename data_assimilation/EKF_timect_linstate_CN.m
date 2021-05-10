function out = EKF_timect_linstate_CN(A,b,g,dAdw,dbdw,dgdx,dgdw,x0,w0,Px0,Pw0,Q,R,yy,uu,dt,opt)

% Extended Kalman Filter - continuous time case (forward euler discretization)
%
% dx/dt = A(w)x + b(w,u) + q      q: N(0,Q)    
% y = g(x,w) + v                  v: N(0,R)   
%
% E[x(0)] = x0, Var[x(0)] = Px0
% E[w] = w0, Var[w] = Pw0
%
% Crank-Nicolson scheme:
% h(x,w,u) = A(w)x + b(w,u)
% (x^{k+1}-x^k)/dt = 1/2 (h(x^{k+1},w,u^{k+1}) + h(x^k,w,u^k))

    opt.dummy = 0;
    
    nX = size(x0,1);
    nW = size(w0,1);

    f_discr = @(x,u1,u2) [(eye(nX)/dt - A(x(nX+1:end))/2)\((eye(nX)/dt + A(x(nX+1:end))/2)*x(1:nX) + (b(x(nX+1:end),u1) + b(x(nX+1:end),u2))/2); ...
                          x(nX+1:end)];

    g_discr = @(x) g(x(1:nX),x(nX+1:end));
    dgdx_discr = @(x) [dgdx(x(1:nX),x(nX+1:end)),dgdw(x(1:nX),x(nX+1:end))];

    x0_discr = [x0;w0];
    P0_discr = [Px0, zeros(nX,nW); zeros(nW,nX), Pw0];
    Q_discr = sqrt(dt)*[Q, zeros(nX,nW); zeros(nW,nX), zeros(nW,nW)]; %TODO: correct normalization of Q and R

    opt.idx_param = nX+1:nX+nW;
    out = EKF_timediscr(f_discr,g_discr,@dfdx_discr,dgdx_discr,x0_discr,P0_discr,Q_discr,R,yy,uu,opt);

    function dfdx = dfdx_discr(X,u1,u2)
        x1 = X(1:nX);
        X2 = f_discr(X,u1,u2);
        x2 = X2(1:nX);
        w = X(nX+1:end);
        
        dxdx = (eye(nX)/dt - A(w)/2)\(eye(nX)/dt + A(w)/2);
        dwdx = zeros(nW,nX);
        dwdw = eye(nW);        
        dxdw = (eye(nX)/dt - A(w)/2)\((dhdw(x2,w,u2)+dhdw(x1,w,u1))/2);
        
        dfdx = [dxdx dxdw; dwdx dwdw];
    end

    function D = dhdw(x,w,u)
        D = dbdw(w,u);
        dAdw_curr = dAdw(w);
        for iW = 1:nW
            D(:,iW) = D(:,iW) + dAdw_curr(:,:,iW)*x;
        end
    end

end