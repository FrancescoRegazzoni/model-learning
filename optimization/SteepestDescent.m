function SteepestDescent(x0,FuncGrad,LinearSearchAlgorithm,PlotFnc,alphaFix,nMaxIter)

    optLS.algorithm = LinearSearchAlgorithm;
    
    it = 0;
    
    PlotFnc(x0);
    
    while (it < nMaxIter)
        TimeIterInit = tic();
        
        it = it+1;
        fprintf('SD iteration %d\n',it)
        
        TimeGradientInit = tic();
        [~,dx] = FuncGrad(x0);
        TimeGradient = toc(TimeGradientInit);
        
        TimeLSInit = tic();
        if LinearSearchAlgorithm == 0
            x0 = x0 - alphaFix*dx;
        else
            [x0,alphaStar] = LinearSearchVectorial(x0,-dx,FuncGrad,FuncGrad,optLS);        
        end
        TimeLS = toc(TimeLSInit);
        
        TimePlotInit = tic();
        x0 = PlotFnc(x0);
        TimePlot = toc(TimePlotInit);
        
        TimeTot = toc(TimeIterInit);
        
        fprintf('   Iter %d   Time: %1.2e   Time Grad: %1.2e   Time LS: %1.2e   Time plot: %1.2e \n', it, TimeTot, TimeGradient, TimeLS, TimePlot)
    end
    
end