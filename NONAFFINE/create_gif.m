function create_gif(filename,pl_state_traj,pl_tau,time_traj,tau,x,pred,Ahat)

desiredFs = round(length(time_traj)/max(time_traj));
for i=1:size(tau,1)
    [tau_resam(i,:), Ty] = resample(tau(i,:),time_traj,desiredFs);
end
x_resam = resample(x,time_traj,desiredFs);

pl_state_traj.XData = Ty(1);
pl_state_traj.YData = x_resam(1);
taucur =  Ahat*tau_resam(:,1);
for i=1:length(pred)
    update_zls(pred(i).pl_zls,pred(i).zls_state_data,...
        pred(i).zls_time_data,pred(i).tr(taucur),pred(i).te(taucur),pred(i).Color)
end

% UPDATE PARAMS PLOTS
tau_hat = tau_resam(:,1);
for i=1:length(tau_hat)
    pl_tau(i).XData = Ty(1);
    pl_tau(i).YData = tau_hat(i,1);
end
drawnow;

for j=1:length(Ty)
    % UPDATE STATE & TRAJECTORY 
    pl_state_traj.XData(end+1) = Ty(j);
    pl_state_traj.YData(end+1) = x_resam(j);
    

    % UPDATE PREDICATE PLOTS
   
    taucur =  Ahat*tau_resam(:,j);
    for i=1:length(pred)
        update_zls(pred(i).pl_zls,pred(i).zls_state_data,...
            pred(i).zls_time_data,pred(i).tr(taucur),pred(i).te(taucur),pred(i).Color)
    end

    % UPDATE PARAMS PLOTS
    tau_hat = Ahat_inv*taucur;
    for i=1:length(tau_hat)
        pl_tau(i).XData(end+1) = Ty(j);
        pl_tau(i).YData(end+1) = tau_hat(i);
    end
    drawnow;
    exportgraphics(gcf,filename,'Append',true);
end


end