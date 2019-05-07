function [alpha, w, b, stat] = norm1svm_dual(X, y, D)
    % alpha = norm1svm_dual(X, y, D)
    
    [d, n] = size(X);
    n1 = sum(y >= 0);
    n2 = sum(y < 0);
    assert(n1 + n2 == n);
 
    A = X(:, y >= 0);
    B = X(:, y < 0);
    
    % Decision variables
    u = sdpvar(n1, 1);
    v = sdpvar(n2, 1);
    
    obj = norm(A * u - B * v, Inf);
    constr = [sum(u) == 1; u >= 0; ...
        sum(v) == 1; v >= 0; ...
        0 <= u; u <= D; ...
        0 <= v; v <= D];
    
    % Solve
    stat = solvesdp(constr, obj);
    u = double(u);
    v = double(v);
    
    % Construct corresponding variables
    alpha = zeros(n, 1);
    alpha(y >= 0) = u;
    alpha(y < 0) = v;
    
    % Solve for primal solution
    [w,b] =  norm1svm_primal_from_dual(X, y, alpha);
    
    %fprintf('Dual objective = %g\n', double(obj));
    stat.obj = double(obj);
end