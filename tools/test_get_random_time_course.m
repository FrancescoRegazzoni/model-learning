% tt = linspace(0,1,100);
tt = 0:1e-2:1;
opt.do_plot = 1;
% opt.time_scale = .1;
opt.time_scale = 5e-2;
opt.dim = 1;
% opt.regDegree = 3;
% opt.fading_coeff = [1 1 .1];
% opt.umin = [-1 3 10];
% opt.umax = [1 4 20];
% opt.method = 'regular_increments';
opt.method = 'coloured_noise';
opt.umin = 0;
opt.umax = 1;
% opt.dudtmin = -8;
% opt.dudtmax = 8;
get_random_time_course(tt,opt);

%axis([0 1 -1 1])