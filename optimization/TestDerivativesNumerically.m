function TestDerivativesNumerically(x0,Evaluation)

ep = 1e-6;

[E,DE] = Evaluation(x0,1,1,0,0,0);
[F,DF] = Evaluation(x0,0,0,1,1,0);
        
DEbis = DF'*F; 

fprintf('Error Gradient (DE - DF''*F): %e\n',norm(DE(:)-DEbis(:),2)/norm(DE(:),2));

fprintf('==== Comparison numerical derivatives (eps = %e) ====\n',ep);

nx = length(x0);
nF = length(F);

DEnum = zeros(nx,1);
DFnum = zeros(nF,nx);
for ix = 1:nx
    incr = zeros(nx,1);
    incr(ix) = ep;
    xnew = x0+incr;    
    Enew = Evaluation(xnew,1,0,0,0,0);
    Fnew = Evaluation(xnew,0,0,1,0,0);
    DEnum(ix) = (Enew-E)/ep;
    DFnum(:,ix) = (Fnew-F)/ep;
end

fprintf('Error Gradient (DE   ): %e\n',norm(DE(:)-DEnum(:),2)/norm(DEnum(:),2));
fprintf('Error Gradient (DF''*F): %e\n',norm(DEbis(:)-DEnum(:),2)/norm(DEnum(:),2));
fprintf('Error Jacobian (DF   ): %e\n',norm(DF(:)-DFnum(:),2)/norm(DFnum(:),2));

subplot(1,2,1)
plot(DE,'linewidth',1.5)
hold on
plot(DEbis)
plot(DEnum,'--')
legend('DE','DEbis','DEnum');

subplot(1,2,2)
plot(DF(:)')
hold on
plot(DFnum(:)')
legend('DF','DFnum');
