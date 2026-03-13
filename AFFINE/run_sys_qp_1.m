%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                               INITIALIZE 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clc; clear; close all;
odeopts = odeset('RelTol',1e-3,'AbsTol',1e-4);
odeopts_pl = odeset('RelTol',1e-5,'AbsTol',1e-5,'Refine',5);
RGB = orderedcolors("gem12");
haxis_lims = [-0.75, 2.25];
filename = 'results_qp_1.gif';
savegif = 0;

% TIME DOMAIN BOUNDS
tmax = 20;

% DYNAMICS
a = 1;
g = @(x) a*(0.5*x + 1.0);
f = @(x) 0.1*tanh(x);
ubound = 0.5;

% TIME STL PARAMS 
t01 = 0;
t02 = 25;

t11 = 3.0;
t12 = 4.0;

t21 = 1.0;
t22 = 2.0;

t31 = 1.0;
t32 = 2.0;

t41 = 8.0;
t42 = 10.0;

t51 = 0;
t52 = 0.5;

title_text = sprintf('$ G_{[%2.1f,%2.1f]} \\left(  F_{[%2.1f,%2.1f]} \\left( \\mu_1 U_{[%2.1f,%2.1f]} \\left( F_{[%2.1f,%2.1f]} \\mu_2 \\right)\\right) \\right)$',...
                        t01,t02,t11,t12, t21,t22, t31,t32);

% PREDICATE FUNCTIONS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pred(1).xc = 1.0;
pred(1).radius = 0.25;
pred(1).h = @(x) 10*(pred(1).radius^2 - (x - pred(1).xc).^2);

pred(1).a0 = t01 + t11;
pred(1).a0_list = [t11; t11];
pred(1).a0_ini = pred(1).a0;
pred(1).a0_list_ini = pred(1).a0_list;

pred(1).tr_list_depends = [1,0,0,0,0,0;
                           1,0,0,0,0,0];
pred(1).a1 = pred(1).tr_list_depends(end,:);
pred(1).tr = @(tau) pred(1).a0 + pred(1).a1*tau;
pred(1).tr_list = @(tau)  pred(1).a0_list + pred(1).tr_list_depends*tau;

pred(1).b0 = t01 + t11 + t21;
pred(1).b0_list = [t11; t11 + t21];
pred(1).b0_ini = pred(1).b0;
pred(1).b0_list_ini = pred(1).b0_list;

pred(1).te_list_depends = [1,0,0,0,0,0;
                           1,1,0,0,0,0];
pred(1).b1 = pred(1).tr_list_depends(end,:);
pred(1).te = @(tau) pred(1).b0 + pred(1).b1*tau;
pred(1).te_list = @(tau) pred(1).b0_list + pred(1).te_list_depends*tau;

pred(1).flag_r = [0;0];
pred(1).flag_e = [0;0];

pred(1).dtr_dtau = [1;0;0;0;0;0];
pred(1).Color = RGB(1,:);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pred(2).xc = 0;
pred(2).radius = 0.25;
pred(2).h = @(x) 10*(pred(2).radius^2 - (x - pred(2).xc).^2);

pred(2).a0 = t01 + t11 + t21 + t31;
pred(2).a0_list = [t11; t11 + t21; t11 + t21 + t31];
pred(2).a0_ini = pred(2).a0;
pred(2).a0_list_ini = pred(2).a0_list;

pred(2).tr_list_depends = [0,0,1,0,0,0;
                           0,0,1,1,0,0;
                           0,0,1,1,1,0];
pred(2).a1 = pred(2).tr_list_depends(end,:);
pred(2).tr = @(tau) pred(2).a0 + pred(2).a1*tau;
pred(2).tr_list = @(tau) pred(2).a0_list + pred(2).tr_list_depends*tau;

pred(2).b0 = t01 + t11 + t21 + t31;
pred(2).b0_list =[t11; t11 + t21; t11 + t21 + t31];
pred(2).b0_ini = pred(2).b0;
pred(2).b0_list_ini = pred(2).b0_list;

pred(2).te_list_depends = [0,0,1,0,0,0;
                           0,0,1,1,0,0;
                           0,0,1,1,1,0];
pred(2).b1 = pred(2).tr_list_depends(end,:);
pred(2).te = @(tau) pred(2).b0 + pred(2).b1*tau;
pred(2).te_list = @(tau) pred(2).b0_list + pred(2).te_list_depends*tau;

pred(2).flag_r = [0;0;0];
pred(2).flag_e = [0;0;0];

pred(2).dtr_dtau = [0;0;1;1;1;0];
pred(2).Color = RGB(2,:);


% ADJACENCY MATRIX
Ahat = [1,0,0,0;
        0,1,0,0;
        1,0,0,0;
        0,1,0,0;
        0,0,1,0;
        0,0,0,1];
Ahat_inv = Ahat'./repmat(sum(Ahat,1)',1,size(Ahat,1));
lb_tau = [0;0;0;0];
ub_tau = [t12 - t11; t22 - t21; t32 - t31; t42 - t41];

tau_hat = ub_tau;
% tau_hat = 0.5*(lb_tau + ub_tau);

tau0 = Ahat*tau_hat;

% PLOT
figure(1); subplot(2,1,1); hold on; grid on; box on;
ylabel('State $x$','Interpreter','latex','FontSize',14);
xlabel('Time $t$','Interpreter','latex','FontSize',14);
hfig = gca(); hfig.TickLabelInterpreter = 'latex';

% RUN ZLS 
legend_text = {};
for i=1:length(pred)
    [pred(i).pl_zls,pred(i).zls_time_data,pred(i).zls_state_data,pred(i).zls] = compute_zls(f,g,pred(i).xc,pred(i).radius,ubound,tmax,odeopts_pl);
    legend_text{end+1} = sprintf("$\\partial V_{h_{%d}}(x,t)$",i);
end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                           SIMULATION SET UP
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% FINAL SIMULATION TIME
tfinal = 30;

% INITIAL CONDITIONS
x0 = -0.6;
t0 = 0;

% CURRENT VALUES
tcur = t0;
xcur = x0;
ucur = 0;
taucur = tau0;

% UPDATE PLOTS
for i=1:length(pred)
    update_zls(pred(i).pl_zls,pred(i).zls_state_data,pred(i).zls_time_data,...
        pred(i).tr(tau0),pred(i).te(tau0),pred(i).Color)
    pl_pred = plot_predicate(pred(i).zls,t0,tfinal,pred(i).Color);
    legend_text{end+1} = sprintf("$h_{%d}\\geq 0$",i);
end
axis([t0,tfinal, haxis_lims(1),haxis_lims(2)])

% PLOT CODE
pl_state_ini = plot(t0,x0,'mo','MarkerSize',8,'MarkerFaceColor','m');
legend_text{end+1} = "Initial State";
pl_state_traj = plot(t0,x0,'m-','LineWidth',3.5);
legend_text{end+1} = "Trajectory";
state_traj = x0;
time_traj = t0;

% SIMULATION STEP SIZE
dt = 1e-2;
% MAX STATE STEP SIZE 
max_dx = 5e-3;
% VARIABLES FOR PLOTTING 
iter = 1;
plotiter = 20;
L1 = legend(legend_text,'Location','northeast','NumColumns',3);
L1.Interpreter = 'latex'; L1.FontSize = 14;
T = title(title_text);
T.Interpreter = 'latex';
T.FontSize = 16;

% PARAMS PLOTS
figure(1); subplot(2,1,2); hold on; grid on; box on;
ylabel('Parameters $\hat{\tau}$[s]','Interpreter','latex','FontSize',14);
xlabel('Time $t$[s]','Interpreter','latex','FontSize',14);
hfig = gca(); hfig.TickLabelInterpreter = 'latex';
legend_text2 = {};
for i=1:length(tau_hat)
    pl_tau(i) = plot(t0,tau_hat(i),'LineStyle','-','LineWidth',2.0,'Color',RGB(4 + i,:));
    legend_text2{end+1} = sprintf("$\\hat{\\tau}_{%d}$",i);
end
axis([t0,tfinal, min(lb_tau), max(ub_tau)]);
L2 = legend(legend_text2,'Location','northeast','NumColumns',length(tau_hat));
L2.Interpreter = 'latex'; L2.FontSize = 16;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                           MAIN SIMULATION LOOP
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
while tcur < tfinal
    
    % COMPUTE CONTROLLER INPUT
    [ucur,omegac] = ucon(xcur,tcur,taucur,Ahat,Ahat_inv,lb_tau,ub_tau,f,g,pred,ubound,odeopts);
    dx = (f(xcur) + g(xcur)*ucur);

    % CONTROL STATE STEP SIZE & UPDATE DT
    if abs(dx*dt)>= max_dx
        dtc = max_dx/abs(dx);
    else
        dtc = dt;
    end
    if max(omegac)*dt>= max_dx
        dtc = max_dx/max(omegac);
    end

    % UPDATE STATE
    xprev = xcur;
    xcur = xcur + dx*dtc;
    taucur = taucur + omegac*dt;

    % CHECK SATISFATION
    for j=1:length(pred)
        [sat_flag,lb_tau,ub_tau,pred(j).flag_r,pred(j).flag_e] = check_sat(tcur,xprev,xcur,pred(j),lb_tau,ub_tau,taucur,Ahat_inv,j);
    end

    % DISPLAY 
    fprintf("Time = %3.2f of %3.2f , Step Size = %3.2e, Prd Vals: %2.3e | %2.3e | %2.3e \n",tcur,tfinal,dtc,pred(1).h(xcur),pred(2).h(xcur),pred(3).h(xcur));


    % CUSTOM FOR REPEAT
    if sat_flag(2)
            beta_curr = pred(2).te(taucur);
            beta_prime = t02;
            if beta_curr<beta_prime
                pred(1).a0 = pred(1).a0_ini + beta_curr;
                pred(1).a0_list = pred(1).a0_list_ini + beta_curr;
              
                pred(1).tr = @(tau) pred(1).a0 + pred(1).a1*tau;
                pred(1).tr_list = @(tau)  pred(1).a0_list + pred(1).tr_list_depends*tau;
                
                pred(1).b0 = pred(1).b0_ini + beta_curr;
                pred(1).b0_list = pred(1).b0_list_ini + beta_curr;

                pred(1).te = @(tau) pred(1).b0 + pred(1).b1*tau;
                pred(1).te_list = @(tau) pred(1).b0_list + pred(1).te_list_depends*tau;
                
                pred(1).flag_r = 0*pred(1).flag_r;
                pred(1).flag_e = 0*pred(1).flag_e;
                % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                pred(2).a0 = pred(2).a0_ini + beta_curr;
                pred(2).a0_list = pred(2).a0_list_ini + beta_curr;
              
                pred(2).tr = @(tau) pred(2).a0 + pred(2).a1*tau;
                pred(2).tr_list = @(tau)  pred(2).a0_list + pred(2).tr_list_depends*tau;
                
                pred(2).b0 = pred(2).b0_ini + beta_curr;
                pred(2).b0_list = pred(2).b0_list_ini + beta_curr;

                pred(2).te = @(tau) pred(2).b0 + pred(2).b1*tau;
                pred(2).te_list = @(tau) pred(2).b0_list + pred(2).te_list_depends*tau;
                
                pred(2).flag_r = 0*pred(2).flag_r;
                pred(2).flag_e = 0*pred(2).flag_e;

                lb_tau = [0;0;0;0];
                ub_tau = [t12 - t11; t22 - t21; t32 - t31; t42 - t41];

            end
    end

    % UPDATE TIME
    tcur = tcur + dtc;

    % SAVE DATA
    state_traj(end+1) = xcur;
    time_traj(end+1) = tcur;

    % UPDATE PLOTS
    if mod(iter,plotiter) == 0|| i==1
        % UPDATE STATE & TRAJECTORY 
        pl_state_traj.XData(end+1) = tcur;
        pl_state_traj.YData(end+1) = xcur;

        % UPDATE PREDICATE PLOTS
        for i=1:length(pred)
            update_zls(pred(i).pl_zls,pred(i).zls_state_data,...
                       pred(i).zls_time_data,pred(i).tr(taucur),pred(i).te(taucur),pred(i).Color)
        end

        % UPDATE PARAMS PLOTS
        tau_hat = Ahat_inv*taucur;
        for i=1:length(tau_hat)
            pl_tau(i).XData(end+1) = tcur;
            pl_tau(i).YData(end+1) = tau_hat(i);
        end
        drawnow;
        if savegif
             exportgraphics(gcf,filename,'Append',true);
        end
        % pause(0.05);
    end
    % INCREMENT ITERATIONS FOR PLOTTING 
    iter=iter+1;

    % COMPUTE FINAL TIME
    for j=1:length(pred)
        tfinal_cands(j) = pred(j).te(taucur);
    end    
    tfinal = max([tfinal_cands,t02]);
end

% UPDATE PLOTS WITH ALL DATA
pl_state_traj.XData = time_traj;
pl_state_traj.YData = state_traj;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                               FUNCTIONS
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% SATISFACTION CHECK
function [sat_flag,lb_tau,ub_tau,flag_r,flag_e] = check_sat(t,xprev,xcur,pred,lb_tau,ub_tau,tau,Ahat_inv,j)
    sat_flag(j) = 0;
    flag_r = pred.flag_r;
    flag_e = pred.flag_e;

    tau_hat = Ahat_inv*tau;

    eps = 1e-2;
    hprev = pred.h(xprev);
    hcur = pred.h(xcur);
    % CHECK PREDICATE REACHING
    if hprev<=0 && hcur>=0
        % CHECK TIMINGS 
        tr = pred.tr_list(tau);
        % RUN THROUGH PATH TIMINGS
        % REACH TIMES
        for i=1:length(tr)
            % IF WITHIN TIME WINDOW
            if abs(t - tr(i))<=eps
                % FIND AND FIX DEPENDENT PARAMETERS
                tau_hat_flag = Ahat_inv*pred.tr_list_depends(i,:)'>0;
                lb_tau(tau_hat_flag) = tau_hat(tau_hat_flag);
                ub_tau(tau_hat_flag) = tau_hat(tau_hat_flag);
                flag_r(i) = 1;
            end
        end  
    end
    % CHECK PREDICATE SATISFACTION
    if hcur>=0
        % CHECK TIMINGS 
        te = pred.te_list(tau);
        % RUN THROUGH PATH TIMINGS
        % END TIMES
        for i=1:length(te)
            % IF WITHIN TIME WINDOW
            if t >= te(i) - eps && ~flag_e(i)
                % FIND AND FIX DEPENDENT PARAMETERS
                tau_hat_flag = Ahat_inv*pred.te_list_depends(i,:)'>0;
                lb_tau(tau_hat_flag) = tau_hat(tau_hat_flag);
                ub_tau(tau_hat_flag) = tau_hat(tau_hat_flag);
                flag_e(i) = 1;
                if i==length(te)
                    sat_flag(j) = 1;
                end
            end
        end  
    end
end


% CONTROLLER 
function [uc,omegac] = ucon(x,t,tau,Ahat,Ahat_inv,lb_tau,ub_tau,f,g,pred,ubound,odeopts)
    num_inps = 1;
    % CONTROLLER PARAMS
    qpopts = optimoptions('quadprog','Display','none','ConstraintTolerance',1e-4);
    clk_1 = 10.0;
    clk_2 = 10.0;
    fact_om = 1.0;
    uref = 0;
    omegaref = -1.0*Ahat_inv*tau;

    % COMPUTE VALUE & GRADIENT
    vf = [];
    dvf_x = [];
    dvf_t = [];
    dvf_dtau = [];
    for i=1:length(pred)
        if t< pred(i).te(tau)
            vf(end+1,1) = value_fun(f,g,pred(i).h,pred(i).xc,pred(i).radius,ubound,odeopts,x,t-pred(i).tr(tau));
            [dvf_x(end+1,1), dvf_t(end+1,1)] = value_fun_grad(f,g,pred(i).h,pred(i).xc,pred(i).radius,ubound,odeopts,x,t-pred(i).tr(tau));
            dvf_dtau(end+1,:) = dvf_t(end,1)*pred(i).dtr_dtau;
        else
            vf(end+1,1) = 0;
            dvf_x(end+1,1) = 0;
            dvf_t(end+1,1) = 0;
            dvf_dtau(end+1,:) = zeros(1,length(tau));
        end
    end

    num_ind = size(Ahat,2);
    lb = [-ubound;-inf(num_ind,1)];
    ub = [+ubound; inf(num_ind,1)];

    H = 2*eye(length(lb));    
    F = -2*[uref; omegaref];

    H(num_inps+1:end,num_inps+1:end) = fact_om*H(num_inps+1:end,num_inps+1:end);
    F(num_inps+1:end) = fact_om*F(num_inps+1:end);

    
    temp_der = -dvf_dtau*Ahat;
    for i=1:length(lb_tau)
        if abs(lb_tau(i) - ub_tau(i))<1e-3
            temp_der(:,i) = zeros(size(temp_der(:,i)));
            lb(num_inps+i) = 0;
            ub(num_inps+i) = 0;
        end
    end

    A = -[dvf_x*g(x), temp_der];
    B = (dvf_x*f(x)+dvf_t)+clk_1.*(vf);

    
    lb_hat = Ahat*lb_tau;
    ub_hat = Ahat*ub_tau;

    Om = -(Ahat_inv*tau - lb_tau).*(Ahat_inv*tau - ub_tau);
    dOm = diag(-(tau - lb_hat) -(tau - ub_hat));

    A = [A; zeros(num_ind,1),-Ahat_inv*dOm*(Ahat)];
    B = [B; clk_2*Om];

    % COMPUTE INPUT
    [sol,~,exitflag] = quadprog(H,F,A,B,[],[],lb,ub,[],qpopts);

    if exitflag==-2
        error("PROBLEM INFEASIBLE")
    else
        uc = sol(1);
        omegac = Ahat*sol(2:end);
    end
end
