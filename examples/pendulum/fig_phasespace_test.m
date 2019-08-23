clear
%% Options
do_save = 1;
print_labels = 0;
showThirdRow = 1;
filename = 'fig/phaseSpaceTest';

lw2 = 1.5;
lw1 = 1.5;

legend_position_bias = [0 -.015 0 0];

%% problem definition
problem = problem_get('pendulum','pendulum.ini');
HFmod = problem.get_model(problem);

%% Test definition
Tmax = 20;
test.tt = [0 Tmax];
test.uu = @(t) [sin(1*t).*cos(1.3*t); ...
                cos(1.8*t).*sin(1*t)];

%% ANN models
load_ANNs;

%% Execution
figure('units','pixel','outerposition',[100,-100,1200,400+200*showThirdRow]);
set(gcf,'Color','w')

nANNmodels = length(ANNmodels_N2);
nCols = 1 + nANNmodels;
nRows = 2 + showThirdRow;
opt_solve.save_x = 1;

outHF = model_solve(test,HFmod,opt_solve);
w = outHF.xx;

subplot(nRows,nCols,1)
plot(outHF.tt,w,'linewidth',lw2)
set(gca,'Color',[1 1 1]*.8);
leg = legend('\theta','d\theta/dt','Orientation','horizontal','Location','south');
leg.Position = leg.Position + legend_position_bias;
legend('boxoff')
axis([0 20 -1 1])  
title('HF model')

subplot(nRows,nCols,nCols+1)
plot(w(1,:),w(2,:),'linewidth',lw1)
if print_labels
    xlabel('\theta')
    ylabel('d\theta/dt')
end
axis([-.3 .3 -1 1])
set(gca,'Color',[1 1 1]*.8);

if showThirdRow
    w = decouple(w);

    subplot(nRows,nCols,2*nCols+1)
    plot(w(1,:),w(2,:),'linewidth',lw1)
    if print_labels
    xlabel('\theta')
    ylabel('d\theta/dt')
    end
        axis([-.3 .3 -3 3])
    set(gca,'Color',[1 1 1]*.8);
end
    
set(gcf, 'InvertHardcopy', 'off')
    
%%
for iNet = 1:nANNmodels
    outANN = model_solve(test,ANNmodels_N2{iNet},opt_solve);
    x = outANN.xx;
    
    subplot(nRows,nCols,1+iNet)
    plot(outANN.tt,x,'linewidth',lw2) 
    leg = legend('x_1','x_2','Orientation','horizontal','Location','south');
    leg.Position = leg.Position + legend_position_bias;
    legend('boxoff') 
    axis([0 20 -1 1])  
    title(sprintf('ANN model #%d',iNet))

    subplot(nRows,nCols,nCols+1+iNet)
    plot(x(1,:),x(2,:),'linewidth',lw1)
    if print_labels
        xlabel('x_1')
        ylabel('x_2')
    end
    axis([-.3 .3 -1 1])
    
    if showThirdRow
        x = decouple(x);        

        % invert sign for visualization purposes
        if iNet == 1 || iNet == 2
            x(2,:) = -x(2,:);
        end

        subplot(nRows,nCols,2*nCols+1+iNet)
        plot(x(1,:),x(2,:),'linewidth',lw1)
        if print_labels
            xlabel('x_1')
            ylabel('x_2')
        end
        axis([-.3 .3 -3 3])
    end

end

if (do_save) 
    create_directory_if_not_found('fig')
    print(filename,'-depsc'); 
end

function x = decouple(x)
    x(2,:) = x(2,:) - mean(x(2,:));
    x(2,:) = x(2,:) - x(1,:) * (x(2,:)*x(1,:)')/(x(1,:)*x(1,:)');
    x(2,:) = x(2,:) / sqrt(mean(x(2,:).^2));
end