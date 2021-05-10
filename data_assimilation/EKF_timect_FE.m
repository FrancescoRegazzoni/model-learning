function out = EKF_timect_FE(f,g,dfdx,dgdx,x0,P0,Q,R,yy,uu,dt,opt)

% Extended Kalman Filter - continuous time case (forward euler discretization)
%
% dx/dt = f(x,u) + w      w = S * <white noise>, Q = SS'
% y = g(x) + v            v = Z * <white noise>, R = ZZ'
%
% E[x(0)] = x0, Var[x(0)] = P0
%
% Dimensions: f [X/T], g [Y], dfdx [1/T], dgdx [Y/X], x0 [X], P0 [X^2]
%             Q [X^2/T], R [Y^2 T]

    opt.dummy = 0;

    f_discr = @(x,u1,u2) x+dt*f(x,u1);
    dfdx_discr = @(x,u1,u2) eye(size(x0,1)) + dt*dfdx(x,u1);
    Q_discr = dt*Q;
    R_discr = R/dt;
    
    out = EKF_timediscr(f_discr,g,dfdx_discr,dgdx,x0,P0,Q_discr,R_discr,yy,uu,opt);
    
end