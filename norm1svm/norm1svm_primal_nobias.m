function w = norm1svm_primal_nobias(X, y, C)
    % w = norm1svm_primal_nobias(X, y, C)

    [d, n] = size(X);
    A = repmat(y', d, 1) .* X;
    
    % Decision variables
    w = sdpvar(d, 1);
    xi = sdpvar(n, 1);
    
    % Build problem
    obj = norm(w,1) + C * sum(xi);
    %obj = norm(w,2) + C * sum(xi);
    constr = [ones(n,1) - xi <= A'*w; xi >= 0];
   
    % Solve
    solvesdp(constr, obj);
    w = double(w);
    
    % Double check hinge-loss errors
    fprintf('err = %g\n', sum(double(xi)))
    terr = zeros(n,1);
    for i=1:n
        terr(i) = min(0, y(i) * (w' * X(:, i)));
    end
    fprintf('terr = %g\n', terr(i));
end