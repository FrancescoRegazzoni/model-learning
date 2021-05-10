function path_out = video_generate(ANNmod, test, opt)
    % Generates the frames to make a video of the model 'ANNmod' for the
    % test case 'test'. 
    % NB: works only for internal dynamics models.

    %% options
    if ~isfield(opt,'do_save')
        opt.do_save = 0;
    end
    if ~isfield(opt,'path_out')
        opt.path_out = ['frames_' datestr(now,'yyyy-mm-dd_hh-MM-ss')];
    end
    if ~isfield(opt,'label_y')
        opt.label_y = 'y';
    end
    if ~isfield(opt,'video_duration')
        opt.video_duration = 1;
    end
    if ~isfield(opt,'fps')
        opt.fps = 25;
    end
    
    %% initialization
    figure('units','pixel','outerposition',[100 100 800 700]);
    margPlotY = .001;
    margPlotX = .01;
    heigths = [.4 .15 .15 .15];
    posY    = [.6 .4  .25 .1];
    
    cmap = get(0, 'DefaultAxesColorOrder');
    col_y_HF    = cmap(1,:);
    col_y_ANN   = cmap(2,:);
    col_u       = cmap(2+1:2+ANNmod.problem.nU,:);
    col_x       = cmap(2+ANNmod.problem.nU+1:2+ANNmod.problem.nU+ANNmod.nX,:);
    
    line_w = 1.2;
    
    lables_HF = {};
    lables_ANN = {};            
    for i = 1:ANNmod.problem.nY
        name_curr = problem_getvariablename(ANNmod.problem,'y',i);
        lables_HF  = [lables_HF  {[name_curr ' (HF)']}];
        lables_ANN = [lables_ANN {[name_curr ' (ANN)']}];
    end
    y_labels = [lables_HF lables_ANN]; 
    
    u_labels = problem_getvariablename_list(ANNmod.problem,'u');
    
    u_min = ANNmod.problem.u_min;
    u_max = ANNmod.problem.u_max;
                
    if opt.do_save
        mkdir(opt.path_out);
    end
        
    dt = ANNmod.dt;
    Tmax = test.tt(end);
    tt = test.tt(1):dt:Tmax;
    nT = length(tt);
    uu = interp_time_series(test.tt,test.uu,tt);
    yy = interp_time_series(test.tt,test.yy,tt);
    
    y_min = min(yy);
    y_max = max(yy);
    
    rawmod = ANNmod.get_raw_model();
    numnF = rawmod.numnF;
    wF = rawmod.wF;
    thetaF = rawmod.thetaF;
    f = @tanh;
    
    period_plot = round(Tmax/opt.video_duration/opt.fps/dt);
    
    % estimation of maximum values for neurons
    maxAlpha = [];
    maxBeta = [];
    maxState = [];
    x0 = ANNmod.x0;
    % time advance
    for iT = 2:round(nT/10)
        input = [uu(:,iT-1);x0];
        [dx,alpha,beta] = EvaluateANN(numnF,wF,thetaF,input,f,1);
        if isempty(maxAlpha)
            maxAlpha = abs(alpha);
        else
            maxAlpha = max(maxAlpha,abs(alpha));
        end
        if isempty(maxBeta)
            maxBeta = abs(beta);
        else
            maxBeta = max(maxBeta,abs(beta));
        end
        if isempty(maxState)
            maxState = abs(x0);
        else
            maxState = max(maxState,abs(x0));
        end
        x0 = x0 + dt*dx;        
    end
    maxState = max(maxState);
    norms = maxBeta;
    norms(1:numnF(1)) = maxAlpha(1:numnF(1));

    axANN = subplot('Position',[0+margPlotX posY(1)+margPlotY 1-2*margPlotX heigths(1)-2*margPlotY]);
    ax1 =   subplot('Position',[0+margPlotX posY(2)+margPlotY 1-2*margPlotX heigths(2)-2*margPlotY]);
    ax2 =   subplot('Position',[0+margPlotX posY(3)+margPlotY 1-2*margPlotX heigths(3)-2*margPlotY]);
    ax3 =   subplot('Position',[0+margPlotX posY(4)+margPlotY 1-2*margPlotX heigths(4)-2*margPlotY]);
    
    %% execution
    idxframe = 1;
    x = zeros(ANNmod.nX,nT);
    x(:,1) = ANNmod.x0;
%     if useG
%         [y(1),~,~] = EvaluateANN(numnG,wG,thetaG,x0,f,BetaOutput);
%     end
    % time advance
    for iT = 2:nT
        input = [uu(:,iT-1);x0];
        [dx,alpha,beta] = EvaluateANN(numnF,wF,thetaF,input,f,1);
        x0 = x0 + dt*dx;
        x(:,iT) = x0;
%         if useG
%             [y(it),~,~] = EvaluateANN(numnG,wG,thetaG,x0,f,BetaOutput);
%         end
        if mod(iT,period_plot) == 0
            subplot(axANN)
            VisualizeANN(numnF,wF,thetaF,alpha,beta,norms);
            title(sprintf('t = %1.2f s',tt(iT)));
            
            subplot(ax1)       
            for i = 1:ANNmod.problem.nU
                u_plots(i) = plot(tt,(uu(i,:)-u_min(i))/(u_max(i)-u_min(i)),'-','linewidth',line_w,'color',col_u(i,:)); hold on;
            end
            axis([0 Tmax -.1 1.1])
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            legend(u_labels,'Orientation','horizontal','location','northeast','AutoUpdate','off')     
            for i = 1:ANNmod.problem.nU
                plot(tt(iT),(uu(i,iT)-u_min(i))/(u_max(i)-u_min(i)),'o','MarkerFaceColor',col_u(i,:),'MarkerEdgeColor','none');
            end
            hold off
                
            subplot(ax2)
            for i = 1:ANNmod.nX
                x_plots(i) = plot(tt(1:iT),x(i,1:iT)/maxState,'linewidth',line_w,'color',col_x(i,:));
                hold on
            end
            axis([0 Tmax -.1 1.1])
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            legendCell = cellstr(num2str((1:ANNmod.nX)', 'x_%d'));
            legend(x_plots, legendCell,'Orientation','horizontal','location','northeast','AutoUpdate','off')
            for i = 1:ANNmod.nX
                plot(tt(iT),x(i,iT)/maxState,'o','MarkerFaceColor',col_x(i,:),'MarkerEdgeColor','none');
            end
            hold off
            
            subplot(ax3)            
            for i = 1:ANNmod.problem.nY
                y_plots(i)                     = plot(tt,yy(i,:),'-','linewidth',line_w,'color',col_y_HF); hold on;
                y_plots(ANNmod.problem.nY + i) = plot(tt(1:iT),x(i,1:iT),'--','linewidth',line_w,'color',col_y_ANN); hold on;
            end
            axis([0 Tmax y_min-.1*(y_max-y_min) y_max+.2*(y_max-y_min)])
            set(gca,'YTickLabel',[]);
            %set(gca,'XTickLabel',[]);
            xlabel('time [s]')
            legend(y_plots, y_labels,'Orientation','horizontal','location','northeast','AutoUpdate','off')
            for i = 1:ANNmod.problem.nY
                plot(tt(iT),yy(i,iT),'o','MarkerFaceColor',col_y_HF,'MarkerEdgeColor','none');
                plot(tt(iT),x(i,iT) ,'o','MarkerFaceColor',col_y_ANN,'MarkerEdgeColor','none');
            end
            hold off
            
            pause(1e-16)
            if opt.do_save
                print(sprintf('%s/%06d',opt.path_out,idxframe),'-dpng');
                idxframe=idxframe+1;
            end
        end
    end
    
    path_out = opt.path_out;
    
end