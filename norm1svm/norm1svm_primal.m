function [w, b, stat] = norm1svm_primal(X, y, C)
    % w = norm1svm_primal(X, y, C)

    [d, n] = size(X);
    A = X * spdiag(y);
    A = [A; y'];
    
    % Decision variables
    w = sdpvar(d, 1);
    xi = sdpvar(n, 1);
    b = sdpvar(1);
    
    % Build problem
    if C == Inf
        obj = norm(w,1);
        constr = ones(n,1) <= A'*[w; b];
    else
        obj = norm(w,1) + C * sum(xi);
        %obj = norm(w,2) + C * sum(xi);
        constr = [ones(n,1) - xi <= A'*[w; b]; xi >= 0];
    end
   
    % Solve
    stat = solvesdp(constr, obj);
    w = double(w);
    b = double(b);
    
    % Double check hinge-loss errors
    %fprintf('err = %g\n', sum(double(xi)))
    terr = zeros(n,1);
    for i=1:n
        terr(i) = max(0, -y(i) * (w' * X(:, i) + b));
    end
    %fprintf('terr = %g (worst = %g)\n', sum(terr), max(terr));
    
    stat.obj = double(obj);
end