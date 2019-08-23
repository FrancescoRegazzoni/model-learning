function res = get_norm_L2_time(tt,yy,idxEval)
    % Computes the L2 norms in time of a time series:
    % \sqrt( \int_0^T |y(t)|^2 dt )
    
    if nargin == 2
        idxEval = 1:length(tt);
    end
    
    ttEval = tt(idxEval);
    dtEval = [ttEval(2:end)-ttEval(1:end-1) 0];
    dtEval(end) = dtEval(end) + (tt(end)-ttEval(end));
    dtEval(1)   = dtEval(1)   + (ttEval(1)-tt(1));
    
    %nY = size(yy,1);
    %T = tt(end)-tt(1);
    
%     res = sqrt(sum(sum((yy(:,idxEval).^2).*dtEval))/nY/T);
%     res = sqrt(sum(sum(yy(:,idxEval).^2,1).*dtEval)/nY/T);
    res = sqrt(sum(sum(yy(:,idxEval).^2,1).*dtEval));
end