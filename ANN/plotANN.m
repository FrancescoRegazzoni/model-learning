function plotANN(numn,w,theta,f,intervals,numpoints)

if numn(1) == 1
    x = linspace(intervals(1),intervals(2),numpoints);
    y = zeros(numn(end),numpoints);
    for ipoint = 1:numpoints
        [output,~,~] = EvaluateANN(numn,w,theta,x(ipoint),f);
        y(:,ipoint) = output;
    end
    plot(x,y);
elseif numn(1) == 2
    x1 = linspace(intervals(1,1),intervals(1,2),numpoints(1));
    x2 = linspace(intervals(2,1),intervals(2,2),numpoints(2));
    [X1,X2] = meshgrid(x1,x2);
    y = zeros(numpoints(1),numpoints(2),numn(end));
    
    for i1point = 1:numpoints(1)        
        for i2point = 1:numpoints(2)
            [output,~,~] = EvaluateANN(numn,w,theta,[x1(i1point);x2(i2point)],f);
            y(i1point,i2point) = output;
        end        
    end
    
    for iOut = 1:numn(end)
        contourf(X1,X2,y(:,:,iOut));
        colorbar()
        pause(.5)
    end
end

