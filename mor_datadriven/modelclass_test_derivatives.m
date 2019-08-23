function modelclass_test_derivatives(modelclass,N,N_alpha,nU)
    x = rand(N,1);
    a = rand(N_alpha,1);
    u = rand(nU,1);
    w = modelclass.paramsF0;
    
    f = modelclass.eval_f(x,u,a,w,1);    
    [Dw,Dx,Da] = modelclass.eval_sensitivity_f();
    
    ep = 1e-6;
    
    for i = 1:N
        x2 = x;
        x2(i) = x(i) + ep;
        f2 = modelclass.eval_f(x2,u,a,w,1);
        Dx_num(:,i) = (f2-f)/ep;
    end
    for i = 1:N_alpha
        a2 = a;
        a2(i) = a(i) + ep;
        f2 = modelclass.eval_f(x,u,a2,w,1);
        Da_num(:,i) = (f2-f)/ep;
    end
    for i = 1:length(w)
        w2 = w;
        w2(i) = w(i) + ep;
        f2 = modelclass.eval_f(x,u,a,w2,1);
        Dw_num(:,i) = (f2-f)/ep;
    end
    
fprintf('==== Comparison numerical derivatives (eps = %e) ====\n',ep);
fprintf('Error Dw: %e\n',norm(Dw-Dw_num,2)/norm(Dw,2));
fprintf('Error Dx: %e\n',norm(Dx-Dx_num,2)/norm(Dx,2));
fprintf('Error Da: %e\n',norm(Da-Da_num,2)/norm(Da,2));

end