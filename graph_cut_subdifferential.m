function drange = graph_cut_subdifferential(W, x)
    % Very naive way to calculate subdifferential
    % Should be able to do this faster with the incidence matrix

    n = length(x);
    
    dmin = zeros(n, 1);
    dmax = zeros(n, 1);
    
    [i, j, v] = find(W);
    E = length(i);
    
    for e=1:E
        % term |xi - xj]
        
        if x(i(e)) == x(j(e))
            dmin(i(e)) = dmin(i(e)) - v(e);
            dmax(i(e)) = dmax(i(e)) + v(e);
            
            dmin(j(e)) = dmin(j(e)) - v(e);
            dmax(j(e)) = dmax(j(e)) + v(e);
        elseif x(i(e)) < x(j(e))
            dmin(i(e)) = dmin(i(e)) - v(e);
            dmax(i(e)) = dmax(i(e)) - v(e);
            
            dmin(j(e)) = dmin(j(e)) + v(e);
            dmax(j(e)) = dmax(j(e)) + v(e);
        else
            dmin(i(e)) = dmin(i(e)) + v(e);
            dmax(i(e)) = dmax(i(e)) + v(e);
            
            dmin(j(e)) = dmin(j(e)) - v(e);
            dmax(j(e)) = dmax(j(e)) - v(e);
        end
    end
    
    drange = [dmin, dmax];
end