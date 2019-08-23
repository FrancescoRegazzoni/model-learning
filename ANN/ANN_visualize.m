function ANN_visualize(numn,w,theta,opt)

if ~isfield(opt,'plot_neurons_vals')
    opt.plot_neurons_vals = 0;
end
if ~isfield(opt,'plot_labels_IO')
    opt.plot_labels_IO = 0;
end

if opt.plot_neurons_vals
    alpha = opt.alpha;
    beta = opt.beta;
    norms = opt.norms;
end

WHratio = 2;

idxW = 1;
A = zeros(sum(numn));
for iLay = 2:length(numn)
    numncumPrec = sum(numn(1:iLay-2));
    numncumCurr = sum(numn(1:iLay-1));
    for iSource=1:numn(iLay-1)
        for iDest=1:numn(iLay)
            A(numncumPrec+iSource,numncumCurr+iDest) = w(idxW);
            idxW = idxW+1;
        end
    end
end

Xpoints = zeros(sum(numn),1);
Ypoints = zeros(sum(numn),1);
iPoint = 1;
for iLay = 1:length(numn)
    for iNodLay = 1:numn(iLay)
        Xpoints(iPoint) = (iLay-1)*WHratio*max(numn)/length(numn);
        Ypoints(iPoint) =  - iNodLay + .5*numn(iLay);        
        iPoint = iPoint+1;
    end
end


if opt.plot_neurons_vals
    A =A+A';
    G = graph(A);
else
    G = digraph(A);
end
    
LWidths = abs(10*G.Edges.Weight)/max(1e-10,max(abs(G.Edges.Weight)));
EdgeCol = zeros(length(w),3);
if opt.plot_neurons_vals
    EdgeCol(:,3) = .9;
else
    EdgeCol(find(G.Edges.Weight>0),2) = .9;
    EdgeCol(find(G.Edges.Weight<0),1) = .9;    
end
%'EdgeLabel',G.Edges.Weight

Labels = cell(sum(numn),1);
if opt.plot_neurons_vals
    Labels(:) = {''};
else
    Labels(1:numn(1)) = {'-'};
    Labels(numn(1)+1:end) = cellstr(num2str(theta));
end
plot(G,'EdgeLabel',{},'LineWidth',LWidths,'NodeLabel',Labels,'XData',Xpoints,'YData',Ypoints,'EdgeColor',EdgeCol)
set(gca, 'visible', 'off');

if opt.plot_neurons_vals
    axis equal
    hold on
    thetvec = linspace(0,2*pi,1e2);
    cosvec = cos(thetvec);
    sinvec = sin(thetvec);
    for iPoint = 1:sum(numn)
        if iPoint <= numn(1)
            val = alpha(iPoint);
            %normVis = normInput(iPoint);
        elseif iPoint > sum(numn)-numn(end)
            val = beta(iPoint);
            %normVis = normOutput(iPoint-(sum(numn)-numn(end)));
        else
            val = beta(iPoint);
            %normVis = normHidden;
        end
        normVis = norms(iPoint);
        R = .5*abs(val)/normVis;
        
        colVal = min(1,max(0,.5 + .5*abs(val)/normVis));
        if val > 0
            col = [.5 colVal .5];
        elseif val < 0
            col = [colVal .5 .5];
        else
            col = [.5 .5 .5];
        end
        if R > 1
           aaaa = 1; 
        end
        patch(Xpoints(iPoint) + R*cosvec,Ypoints(iPoint) + R*sinvec,col,'EdgeColor','none')
    end
    hold off
end

xGap = .5;

if opt.plot_labels_IO
    for i = 1:numn(1)
        iPoint = i;
        text(Xpoints(iPoint)-xGap,Ypoints(iPoint),...
            sprintf('$%s \\rightarrow$',opt.label_in{i}),'Interpreter','latex','HorizontalAlignment', 'right')
    end
    for i = 1:numn(end)
        iPoint = sum(numn(1:end-1))+i;
        text(Xpoints(iPoint)+xGap,Ypoints(iPoint),...
            sprintf('$\\rightarrow %s$',opt.label_out{i}),'Interpreter','latex')
    end
end
