function dataset = get_dataset(name)

    t1 = 1;
    t2 = 13;
    t3 = 25;      
    iTests = 1;
    
    switch name
        case 'train'
            u1_steps = linspace(-1,1,3);
            u2_steps = linspace(-1,1,4);
            for u1_val = u1_steps
                for u2_val = u2_steps    
                    u1 = @(t) u1_val*(t>t1).*(t<t2);
                    u2 = @(t) u2_val*(t>t1).*(t<t2);
                    dataset{iTests}.tt = [0 t3];
                    dataset{iTests}.uu = @(t) [u1(t);u2(t)];
                    iTests = iTests+1;
                end
            end

            u1 = @(t) sin(2*t).*min(1,exp(-(t-8)/6));
            u2 = @(t) cos(1*t).*min(1,exp(-(t-10)/4));
            dataset{iTests}.tt = [0 45];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;

            u1 = @(t) .5*sin(1*t).*min(1,exp(-(t-10)/4));
            u2 = @(t) cos(2*t).*min(1,exp(-(t-10)/5));
            dataset{iTests}.tt = [0 45];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;

            u1 = @(t) (-.4 + .4*cos(.7*t)).*(1 + sin(3*t).^2)/(2*.8).*min(1,exp(-(t-20)/4));
            u2 = @(t) (.1 + .6*cos(.6*t)).*(1 + sin(3.2*t).^2)/(2*.7).*min(1,exp(-(t-19)/3));
            dataset{iTests}.tt = [0 45];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;

            u1 = @(t) sin(2*t).*min(1,exp(-(t-10)/4));
            u2 = @(t) (-.3 + .6*cos(.9*t)).*min(1,exp(-(t-10)/6));
            dataset{iTests}.tt = [0 45];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];

        case 'validation'       
            Tmax = 50;
            u1 = @(t) (.5 + .4*sin(1*t)) .* (t < 5) ...
                + (-.3 + .6*cos(2*t)) .* (t >= 5) .* (t < 16) ...
                + .8 * (t >= 16) .* (t < 32) ...
                - sin(1*t).*exp(-t/10) .* (t >= 32) .* (t < 50);       
            u2 = @(t) (.4 + .6*sin(2*t)) .* (t < 6) ...
                + (-.5 + .3*sin(.8*t)) .* (t >= 6) .* (t < 18) ...
                - ( .2 + .8*sin(1*t).*exp(-t/10) ) .* (t >= 18) .* (t < 30) ...
                - .8 * (t >= 30) .* (t < 50);  
            dataset{iTests}.tt = [0 Tmax];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;

            u1 = @(t) .5*(t>t1+3).*(t<t2+5);
            u2 = @(t) .5*(t>t1).*(t<t2); 
            dataset{iTests}.tt = [0 t3];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;

            u1 = @(t) (-.3 + .5*sin(.6*t)).*(1 + cos(3*t).^2)/(2*.8);
            u2 = @(t) (.2 + .8*cos(.9*t)).*(1 + sin(4*t).^2)/(2*1); 
            dataset{iTests}.tt = [0 45];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;

            u1 = @(t) (.2 + .4*sin((1+.4*sin(.1*t)).*t)).*(1 + cos(.5*t)) + .1;
            u2 = @(t) (-.1 + .6*cos(.3*min(10,exp(t/10)).*t)).*(.5 + .2*sin(.3*t).^2); 
            dataset{iTests}.tt = [0 45];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;

            u1 = @(t) .3+.4*min(1,exp(-(t-10)/3))-.6*min(1,exp(-(t-26)/1))+.4*min(1,exp(-(t-34)/.5))-.3.*min(1,exp(-(t-40)/2));
            u2 = @(t) (.5+.3*sin(.5*t)).*(-.2+.4*min(1,exp(-(t-5)/.7))-.6*min(1,exp(-(t-20)/1.6))+.3*min(1,exp(-(t-32)/.5))-.1*min(1,exp(-(t-42)/1))); 
            dataset{iTests}.tt = [0 45];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;
            
        case 'test'
            t1 = 1;
            t2 = 6;
            t3 = 13; 
            u1_steps = linspace(-1,1,8);
            u2_steps = linspace(-1,1,8);
            for u1_val = u1_steps
                for u2_val = u2_steps    
                    u1 = @(t) u1_val*(t>t1).*(t<t2);
                    u2 = @(t) u2_val*(t>t1).*(t<t2);
                    dataset{iTests}.tt = [0 t3];
                    dataset{iTests}.uu = @(t) [u1(t);u2(t)];
                    iTests = iTests+1;
                end
            end

            u1base = [0 .2 -.3];
            u2base = [.1 .4 -.1];
            u1ampl = [.8 .4 .6];
            u2ampl = [.3 .5 .2];
            freq1 = linspace(.5,3,4);
            freq2 = linspace(.5,3,4);
            for i = 1:length(u1base)
                for i1 = 1:length(freq1)
                    for i2 = 1:length(freq2)
                        u1 = @(t) (u1base(i) + u1ampl(i)*sin(freq1(i1)*t));
                        u2 = @(t) (u2base(i) + u2ampl(i)*sin(freq1(i2)*t - .3)); % so that they are never in phase
                        dataset{iTests}.tt = [0 10];
                        dataset{iTests}.uu = @(t) [u1(t);u2(t)];
                        iTests = iTests+1;
                    end
                end
            end

            Tmax = 50;
            u1 = @(t) (.5 + .4*sin(1*t)) .* (t < 5) ...
                + (-.3 + .6*cos(2*t)) .* (t >= 5) .* (t < 16) ...
                + .8 * (t >= 16) .* (t < 32) ...
                - sin(1*t).*exp(-t/10) .* (t >= 32) .* (t < 50);       
            u2 = @(t) (.4 + .6*sin(2*t)) .* (t < 6) ...
                + (-.5 + .3*sin(.8*t)) .* (t >= 6) .* (t < 18) ...
                - ( .2 + .8*sin(1*t).*exp(-t/10) ) .* (t >= 18) .* (t < 30) ...
                - .8 * (t >= 30) .* (t < 50);  
            dataset{iTests}.tt = [0 Tmax];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;

            u1 = @(t) (-.3 + .5*sin(.6*t)).*(1 + cos(3*t).^2);
            u2 = @(t) (.2 + .8*cos(.9*t)).*(1 + sin(4*t).^2);
            dataset{iTests}.tt = [0 45];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;

            u1 = @(t) (.2 + .4*sin((1+.4*sin(.1*t)).*t)).*(1 + cos(.5*t)) + .1;
            u2 = @(t) (-.1 + .6*cos(.3*min(10,exp(t/10)).*t)).*(.5 + .2*sin(.3*t).^2);
            dataset{iTests}.tt = [0 45];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;

            u1 = @(t) .3+.4*min(1,exp(-(t-10)/3))-.6*min(1,exp(-(t-26)/1))+.4*min(1,exp(-(t-34)/.5))-.3*min(1,exp(-(t-40)/2));
            u2 = @(t) (.5+.3*sin(.5*t)).*(-.2+.4*min(1,exp(-(t-5)/.7))-.6*min(1,exp(-(t-20)/1.6))+.3*min(1,exp(-(t-32)/.5))-.1*min(1,exp(-(t-42)/1)));
            dataset{iTests}.tt = [0 45];
            dataset{iTests}.uu = @(t) [u1(t);u2(t)];
            iTests = iTests+1;
            
            nRand = 10;
            Tmax = 10;
            dt_generation = 1e-2;
            tt_rand = 0:dt_generation:Tmax;
            opt_rnd.dim = 2;
            opt_rnd.umin = -.8;
            opt_rnd.umax = .8;
            opt_rnd.regDegree = 2;
            opt_rnd.time_scale = 1e-1;
            opt_rnd.do_plot = 1;
            rng(0) % Inizialize the random number generator in a predictable way
            for i = 1:nRand
                dataset{iTests}.tt = tt_rand;
                dataset{iTests}.uu = get_random_time_course(tt_rand,opt_rnd);
                iTests = iTests+1;
            end

        otherwise
            error('dataset name not recognized')
    end
end