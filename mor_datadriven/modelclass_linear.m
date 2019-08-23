function modelclass = modelclass_linear(optionsfile,problem,N,N_alpha,useG)
    
    nU = problem.nU;
    nY = problem.nY;
    
    modelclass.nwF = N*(N+nU+N_alpha+1);
    if useG
        modelclass.nwG = nY*(N+1);
    end

    %% Labels
    modelclass.name_modelclass = '_linear';    
    
    %% Renormalization
    modelclass.renormalize_allows = 0;
    
    %% Plot options
    modelclass.plot_x = 1;
    modelclass.x_idx_plot = nY+1:N;
    
    %% Model definition
    curr_Fmatrix = [];
    curr_Gmatrix = [];
    curr_Finput = [];
    curr_Ginput = [];
    
    modelclass.eval_f = @eval_f; 
    function Foutput = eval_f(x,u,a,w,dt)
        curr_Fmatrix = reshape(w,N,N+nU+N_alpha+1);
        curr_Finput = [u;x;a;1];
        Foutput = curr_Fmatrix*curr_Finput;
%         xnew = x + dt*curr_Fmatrix*curr_Finput;
    end

    modelclass.eval_g = @eval_g; 
    function outputG = eval_g(x,w)
        curr_Gmatrix = reshape(w,nY,N+1);
        curr_Ginput = [x;1];
        outputG = curr_Gmatrix*curr_Ginput;
    end

    modelclass.eval_sensitivity_f = @eval_sensitivity_f;
    function [Dw,Dx,Da] = eval_sensitivity_f()
        Dx = curr_Fmatrix(:,nU+1:nU+N);
        Da = curr_Fmatrix(:,nU+N+1:nU+N+N_alpha);
        Dw = diagonal_sensitivity(N,curr_Finput);
    end

    modelclass.eval_sensitivity_g = @eval_sensitivity_g;
    function [Dw,Dx] = eval_sensitivity_g()   
        Dx = curr_Gmatrix(:,1:N);
        Dw = diagonal_sensitivity(nY,curr_Ginput);
    end
      
    modelclass.eval_sensitivity_dfdaPen = @() error('not (yet) implemented');
    modelclass.eval_sensitivity_origStabPen = @() error('not (yet) implemented');

    function S = diagonal_sensitivity(n,x)
        S = [];
        for i = 1:length(x)
            S = [S eye(n)*x(i)];
        end
    end
end