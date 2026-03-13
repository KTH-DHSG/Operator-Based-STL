function [pl_zls,zls_time_data,zls_state_data,zls] = compute_zls(f,g,root_g_sq,xc,radius,ubound,tmax,odeopts_pl)
    % ZERO LEVEL SET 
    zls = [xc + radius, xc - radius];
    xzls = {};
    tzls = {};
    for i=1:length(zls)
        % SOLVE ODE
        [tzlstmp,xzlstmp] = ode89(@(t,x) fmax(f,g,root_g_sq,x,xc,radius,ubound),[0,tmax],zls(i),odeopts_pl);
        xzls{i} = xzlstmp';
        tzls{i} = -tzlstmp';
    end
    
    % PLOT ZLS
    zls_time_data  = tzls;
    zls_state_data = xzls;
    pl_zls = plot([tzls{1},nan,tzls{2}],[xzls{1},nan,xzls{2}],'r-','LineWidth',2.5);
    drawnow;
end