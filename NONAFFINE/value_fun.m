% COMPUTE VALUE FUNCTION POINT-WISE
function vf = value_fun(f,g,root_g_sq,h,xc,radius,ubound,odeopts,x,t)
% SOLVE ODE
if t<0 && abs(x-xc)>radius
    [ttmp,xtmp] = ode45(@(t,y) fmax(f,g,root_g_sq,y,xc,radius,ubound),[0,t],x,odeopts);
    vf = max(h(xtmp));
elseif t>=0
    vf = h(x);
else
    vf = h(xc);
end
end