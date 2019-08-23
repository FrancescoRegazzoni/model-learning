function pen_handler = penalization_params_direct(optionsfile)

    fPen_s              = iniread(optionsfile,'Penalizations','pen_f','s','0');
    fPen_type           = iniread(optionsfile,'Penalizations','pen_f_type','s','q');
    fPen_interval       = iniread(optionsfile,'Penalizations','pen_f_interval','d',[-1,1]);
    gPen_s              = iniread(optionsfile,'Penalizations','pen_g','s','0');
    gPen_type           = iniread(optionsfile,'Penalizations','pen_g_type','s','q');
    gPen_interval       = iniread(optionsfile,'Penalizations','pen_g_interval','d',[-1,1]);
    aPen_s              = iniread(optionsfile,'Penalizations','pen_a','s','0');
    aPen_type           = iniread(optionsfile,'Penalizations','pen_a_type','s','q');
    aPen_interval       = iniread(optionsfile,'Penalizations','pen_a_interval','d',[-1,1]);
    
    pen_handler.num_pen = 3;
    pen_handler.names = {'ERRf','ERRg','ERRa'};
    pen_handler.coefficient_string = {fPen_s, gPen_s, aPen_s};
    pen_handler.dependence_f = [1,0,0];
    pen_handler.dependence_g = [0,1,0];
    pen_handler.dependence_a = [0,0,1];
    pen_handler.dependence_i = [0,0,0];
       
    fPen_func = []; fPen_dfunc = [];
    gPen_func = []; gPen_dfunc = [];
    aPen_func = []; aPen_dfunc = [];
    pen_handler.initialize = @initialize;
    function pen_handler = initialize(pen_handler,problem,N,N_alpha,useG,nS,nwF,nwG,nAL,nIC,misc)
        if pen_handler.pen_active(1)
            [fPen_func,fPen_dfunc] = get_pen_functions(fPen_type,fPen_interval);
            pen_handler.norm_pen(1) = 1/nwF;
            pen_handler.size_pen(1) = nwF;
        end
        if pen_handler.pen_active(2)
            [gPen_func,gPen_dfunc] = get_pen_functions(gPen_type,gPen_interval);
            pen_handler.norm_pen(2) = 1/nwG;
            pen_handler.size_pen(2) = nwG;
        end
        if pen_handler.pen_active(3)
            [aPen_func,aPen_dfunc] = get_pen_functions(aPen_type,aPen_interval);
            pen_handler.norm_pen(3) = 1/nAL;
            pen_handler.size_pen(3) = nAL;
        end
    end

    pen_handler.evaluate = @evaluate;
    function out = evaluate(paramsF,paramsG,Alpha,IC,compRes,compGrad,compJac,pen_handler,numIter,modelclass)
        out.dummy = 0;
        if pen_handler.coef(1)>0, out.errors(1) = sum(fPen_func(paramsF).^2); end
        if pen_handler.coef(2)>0, out.errors(2) = sum(gPen_func(paramsG).^2); end
        if pen_handler.coef(3)>0, out.errors(3) = sum(aPen_func(Alpha).^2); end
        
        if compGrad
            if pen_handler.coef(1)>0, out.gradF{1} = 2*fPen_func(paramsF).*fPen_dfunc(paramsF); end
            if pen_handler.coef(2)>0, out.gradG{2} = 2*gPen_func(paramsG).*gPen_dfunc(paramsG); end
            if pen_handler.coef(3)>0, out.gradA{3} = 2*aPen_func(reshape(Alpha,[],1)).*aPen_dfunc(reshape(Alpha,[],1)); end
        end
        
        if compRes
            if pen_handler.coef(1)>0, out.res{1} = fPen_func(paramsF); end
            if pen_handler.coef(2)>0, out.res{2} = gPen_func(paramsG); end
            if pen_handler.coef(3)>0, out.res{3} = aPen_func(reshape(Alpha,[],1)); end
        end
        
        if compJac
            if pen_handler.coef(1)>0, out.jacF{1} = diag(fPen_dfunc(paramsF)); end
            if pen_handler.coef(2)>0, out.jacG{2} = diag(gPen_dfunc(paramsG)); end
            if pen_handler.coef(3)>0, out.jacA{3} = diag(aPen_dfunc(reshape(Alpha,[],1))); end
        end            
    end

    function [Pen_func,Pen_dfunc] = get_pen_functions(Pen_type,Pen_interval)
        switch Pen_type
            case 'q'
        %             fPen_func = @(x) .5*x.^2;
        %             fPen_dfunc = @(x) x;
                Pen_func = @(x) x;
                Pen_dfunc = @(x) 0*x + 1;
            case 'i'
        %             fPen_func = @(x) .5*((x<fPen_interval(1)).*(x-fPen_interval(1)).^2 + (x>fPen_interval(2)).*(x-fPen_interval(2)).^2);
        %             fPen_dfunc = @(x) (x<fPen_interval(1)).*(x-fPen_interval(1)) + (x>fPen_interval(2)).*(x-fPen_interval(2));  
                Pen_func = @(x) (x<Pen_interval(1)).*(x-Pen_interval(1)) + (x>Pen_interval(2)).*(x-Pen_interval(2));  
                Pen_dfunc = @(x) (x<Pen_interval(1)) + (x>Pen_interval(2));  
        end
    end
end