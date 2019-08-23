tt_in =  -1:.1:4;
tt_out = 0:.5:3;
vv_in = [sin(tt_in);cos(tt_in)];

vv_out_1 = interp_time_series(tt_in,vv_in,tt_out,struct('mode','pointwise'));
vv_out_2 = interp_time_series(tt_in,vv_in,tt_out,struct('mode','mean_forward'));

figure()
plot(tt_in,vv_in,'o-');
hold on
plot(tt_out,vv_out_1,'+-');
plot(tt_out,vv_out_2,'+-');