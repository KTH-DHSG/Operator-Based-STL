%% FUNCTIONS
function dx = fmax(f,g,root_g_sq,x,xc,radius,ubound)

    ucand = [-ubound,ubound];
    if root_g_sq(x)>0
        uprime = sqrt(root_g_sq(x));
        if uprime>= -ubound && uprime <= ubound
            ucand(end+1) = uprime; 
        end
        if -uprime>= -ubound && -uprime <= ubound
            ucand(end+1) = -uprime; 
        end
    end


    gcand = g(x,ucand)*(xc - x);
    [~,id] = max(gcand);
    dx = -(f(x) + g(x,ucand(id)));
    % ucand(id)
end
