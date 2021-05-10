function out = mor_ANN_interpolate_test(problem,in,dt_integration,dt_evaluation,interpolation_mode_u,interpolation_mode_y)
    
    %% getting dimensions
    if isfield(in,'uu')
        nU = size(in.uu,1);
    else
        nU = 0;
    end
    nY = size(in.yy,1);
    
    
    %% setting missing labels
    out.eval_diff = 1;
    if isfield(in,'eval_diff')
        if ~in.eval_diff
            out.eval_diff = 0;
        end
    end
    
    out.closed_loop = 0;
    if isfield(in,'closed_loop')
        if in.closed_loop
            out.closed_loop = 1;
        end
    end
    
    %% time step computation
    out.Tmax = max(in.tt);
    
    if isfield(in,'dt_eval')
        dt_evaluation_curr = in.dt_eval;
    else
        dt_evaluation_curr = dt_evaluation;
    end

    if dt_integration(1) == '@'
        dt_intFunc = dt_integration(2:end);
        Newtt(1) = 0;
        t = 0;
        goon = 1;
        while t<out.Tmax
            dt_curr = eval(dt_intFunc);
            t = min(t+dt_curr,out.Tmax);        
            Newtt = [Newtt t];
        end

        if dt_evaluation_curr(1) ~= 'x'
            error('dt_integration wrong format')
        end
    else
        dt_integrationVal = str2double(dt_integration);
        Newtt = 0:dt_integrationVal:out.Tmax;
        if dt_evaluation_curr(1) ~= 'x'
            dt_evaluationVal = str2double(dt_evaluation_curr);
            perEval = dt_evaluationVal/dt_integrationVal;
            out.idxEval = 1:perEval:size(Newtt,2)-1;
        end
    end

    if dt_evaluation_curr(1) == 'x'
        perEval = str2num(dt_evaluation_curr(2:end));
        out.idxEval = 1:perEval:size(Newtt,2)-1;
    end
    
    if isfield(in,'uu')
        out.uu = interp_time_series(in.tt,in.uu,Newtt,struct('mode',interpolation_mode_u));
    end
    if ~isfield(in,'tt_y')
        in.tt_y = in.tt;
    end
    out.yy = interp_time_series(in.tt_y,in.yy,Newtt,struct('mode',interpolation_mode_y));
    if isfield(in,'yy_ex')
        out.yy_ex = interp_time_series(in.tt_y,in.yy_ex,Newtt,struct('mode',interpolation_mode_y));
    end
        
%     if ~isequal(in.tt,Newtt)        
%         for iU = 1:nU
%             out.uu(iU,:) = interp1(in.tt,in.uu(iU,:),Newtt);
%         end
%         for iY = 1:nY
%             out.yy(iY,:) = interp1(in.tt,in.yy(iY,:),Newtt);
%         end
%         if isfield(in,'yy_ex')
%             for iY = 1:nY
%                 out.yy_ex(iY,:) = interp1(in.tt,in.yy_ex(iY,:),Newtt);
%             end
%         end
%     else
%         if nU > 0
%             out.uu = in.uu;
%         end
%         out.yy = in.yy;
%         if isfield(in,'yy_ex')
%             out.yy_ex = in.yy_ex;
%         end
%     end
    out.tt = Newtt;
    out.nT = size(out.tt,2);
    out.nTeval = length(out.idxEval);
    
    out.dt = out.tt(2:end)-out.tt(1:end-1);
    out.dtEval = [out.tt(out.idxEval(2:end))-out.tt(out.idxEval(1:end-1)) (out.tt(end)-out.tt(out.idxEval(end)))];
    
    if problem.samples_variability && isfield(in, 'alpha')
        out.alpha = in.alpha;
    end
end