ddt = 0.05;
ddx = 0.025;
tmin = -5.0;
xmin = -5.0;
xmax = 10;

figure(); hold on; box on; 
[xg,tg] = meshgrid(xmin:ddx:xmax,tmin:ddt:0);

for i=1:1
    % subplot(1,length(pred),i); hold on; box on; 
    for i1 = 1:size(xg,1)
        for i2 = 1:size(xg,2)
            fprintf("I1: %d, I2: %d \n",i1,i2)
            vf{i}(i1,i2) = value_fun(f,g,root_g_sq,pred(i).h,pred(i).xc,pred(i).radius,ubound,odeopts,xg(i1,i2),tg(i1,i2));

        end
    end
    [pred(i).pl_zls,pred(i).zls_time_data,pred(i).zls_state_data,pred(i).zls] = compute_zls(f,g,root_g_sq,pred(i).xc,pred(i).radius,ubound,tmax,odeopts_pl);
    
    % maxv = max(max(vf{i}));
    % minv = min(min(vf{i}));
    % vf{i} = (vf{i} - maxv);

    s{i} = surf(tg,xg,tanh(0.5*vf{i}));
    s{i}.EdgeAlpha = 0.0;
    shading interp;
    ylabel('State $x$','Interpreter','latex','FontSize',14);
    xlabel('Time $t$','Interpreter','latex','FontSize',14);
    hfig = gca(); hfig.TickLabelInterpreter = 'latex';
    axis tight;
    c =  colorbar;
    c.TickLabelInterpreter = 'latex';
    axis([tmin 0 xmin xmax min(min(vf{i})) max(max(vf{i}))]);

end

