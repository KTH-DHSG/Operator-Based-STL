%% FUNCTIONS
function dx = fmax(f,g,x,xc,radius,ubound)

    ucand = [-ubound,ubound];
    gcand = g(x)*ucand*(xc - x);
    [~,id] = max(gcand);
    dx = -(f(x) + g(x)*ucand(id));
    % ucand(id)
end
