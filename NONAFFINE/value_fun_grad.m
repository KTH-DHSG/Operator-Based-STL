% COMPUTE VALUE FUNCTION POINT-WISE
function [dvf_x, dvf_t] = value_fun_grad(f,g,root_g_sq,h,xc,radius,ubound,odeopts,x,t)
eps = 1e-5;
% SOLVE ODE
% [~,xtmp] = ode45(@(t,y) -fmax(f,g,y,xc,ubound),[0,abs(t)],x,odeopts);
% [~,xtmp_pl_t] = ode45(@(t,y) -fmax(f,g,y,xc,ubound),[0,abs(t) - eps],x,odeopts);
% [~,xtmp_pl_x] = ode45(@(t,y) -fmax(f,g,y,xc,ubound),[0,abs(t)],x + eps,odeopts);
% vf = h(xtmp(end));
% vf_pl_t = h(xtmp_pl_t(end));
% vf_pl_x = h(xtmp_pl_x(end));

vf = value_fun(f,g,root_g_sq,h,xc,radius,ubound,odeopts,x,t);
vf_pl_t = value_fun(f,g,root_g_sq,h,xc,radius,ubound,odeopts,x,t+eps);
vf_pl_x = value_fun(f,g,root_g_sq,h,xc,radius,ubound,odeopts,x+eps,t);

dvf_x = (vf_pl_x - vf)/eps;
dvf_t = (vf_pl_t - vf)/eps;
end