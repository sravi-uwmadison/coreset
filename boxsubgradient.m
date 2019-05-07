function [x] = boxsubgradient(A, b, df, x0, gamma, kmax)
    % x = boxsubgradient(A, b, df, x0)
    %
    %   Solves a nonsmooth problem using Frank-Wolfe, where the subgradient
    %   is always a "box."
    %
    %     min f(x)
    %      st A x <= b
    %
    % df: Subdifferential closure of the form [dmin, dmax] = df(x).
    %     Any subgradient g of f at x will obey dmin <= g <= dmax at all
    %     coordinates.
    %
    % Runs for kmax iterations, or until the primal-dual gap is less than
    % gamma, whichever occurs first
    %
    
    n = length(x0);
    x = x0;
    Am = size(A, 1);
    obj = zeros(kmax,1);    
    for k = 0:kmax-1
        alpha = 2/(k+2);
        
        % Calculate subdifferential.
        drange = df(x);
        dmin = drange(:, 1);
        dmax = drange(:, 2);
        
        % Get search direction
        cplex = Cplex('search-direction');
        
        % z cols first
        cplex.addCols(zeros(n, 1), [], -Inf(n, 1), Inf(n, 1));
        
        % mu cols
        cplex.addCols(ones(n, 1), [], -Inf(n, 1), Inf(n, 1));
        
        % Add in A z <= b
        cplex.addRows(-Inf(Am, 1), [A, sparse(Am, n)], b);
        
        Dmin = spdiag(dmin);
        Dmax = spdiag(dmax);
        In = speye(n);
        cplex.addRows(-Dmin * x, [-Dmin, In], Inf(n, 1));
        cplex.addRows(-Dmax * x, [-Dmax, In], Inf(n, 1));
        
        % Solve search direction problem
        cplex.solve();
        s = cplex.Solution.x(1:n);
        
        x = x + alpha * (s - x);        
        
        
        % TODO: calculate primal-dual gap
    end
end