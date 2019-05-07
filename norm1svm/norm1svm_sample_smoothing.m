function [alpha, w, b, stat] = norm1svm_sample_smoothing(X, y, D, opts)
    % [alpha, w, b, stat] = norm1svm_frank_wolfe(X, y, D, opts)
    
    tic;
    
    
    % Set up options
    defaults.max_iters = Inf;
    defaults.converge_tol = 1e-8;
    defaults.line_search = 1;
    defaults.print_interval = 10;
    opts = optdefaults(opts, defaults);
    
    % Structure for saving optimization progress/status
    stat.objective = [];
    stat.gap = [];
    stat.active_indices = [];
    
    % Construct problem matrices
    [d, n] = size(X);
    n1 = sum(y >= 0);
    n2 = sum(y < 0);
    assert(n1 + n2 == n);
    
    A = X(:, y >= 0);
    B = X(:, y < 0);
    
    [bnd, M, ~] = converge_bounds(X, y);
    lan_bnd_sum = 0;
    
    % Initial point
    u = zeros(n1, 1);
    v = zeros(n2, 1);
    K = ceil(1/D);
    assert(K < n1 && K < n2);
    u(1:K) = 1/K;
    v(1:K) = 1/K;
    
    diamX = 2*sqrt(2); % TODO?
    
    assert(min(u) >= -1e-6 && max(u) <= D + 1e-6 && abs(sum(u) - 1) < 1e-6);
    assert(min(v) >= -1e-6 && max(v) <= D + 1e-6 && abs(sum(v) - 1) < 1e-6);
    
    yalopts = sdpsettings('verbose', 0);
    
    obj = norm(A*u - B*v, Inf);
    fprintf(1, '0 %g\n', obj);
    
    dualobj = -Inf;
    k = 1;
    
    while obj - dualobj > opts.converge_tol && k < opts.max_iters
        if norm(A*u - B*v, Inf) < 1e-8
            break
        end
        
        % Sample the gradient
        grad = zeros(n, 1);
        x = zeros(n, 1);
        x(y >= 0) = u;
        x(y < 0) = v;
        
        Tk = k;
        muk = (d ^ (0.25) * diamX) / sqrt(k);
        for sampi=1:Tk
            % Sample from a ball around the current iterate
            x_xi = x + sample_ball(n, muk);
            
            % Calculate gradient at x_xi
            Ax_diff = X * (y .* x_xi);
            infnrm = norm(Ax_diff, Inf);
            Ax_grad = zeros(d, 1);
            Ax_grad(Ax_diff == infnrm) = 1;
            Ax_grad(Ax_diff == -infnrm) = -1;
            grad_xi = y .* (X' * Ax_grad);
            
            grad = grad + grad_xi;
        end
        grad = grad / Tk;
        
        gu = grad(y >= 0);
        gv = grad(y < 0);
        
        % Set up Frank-Wolfe subproblem for this iteration
        zu = sdpvar(n1, 1);
        zv = sdpvar(n2, 1);
        
        
        constr = [sum(zu) == 1; ...
                  sum(zv) == 1; ...
                  0 <= zu; zu <= D; ...
                  0 <= zv; zv <= D];
        solvesdp(constr, gu' * zu + gv' * zv, yalopts);
        
        zu = double(zu);
        zv = double(zv);
        
        
        assert(min(zu) >= -1e-6 && max(zu) <= D + 1e-6 && abs(sum(zu) - 1) < 1e-6);
        assert(min(zv) >= -1e-6 && max(zv) <= D + 1e-6 && abs(sum(zv) - 1) < 1e-6);
    
        Az = A*zu - B*zv;
        %dualobj = norm(Az, Inf);
        dualobj = gu' * zu + gv' * zv;

        %if dualobj > obj + opts.converge_tol
        %    break
        %end

        z = zeros(n,1);
        z(y >= 0) = zu;
        z(y < 0) = zv;
        
        % TODO: should be able to do exact line search, might improve things
        if opts.line_search
            linesearch = @(x) norm(X * (y .* x), Inf);
            a = binary_line_search(x, z - x, linesearch, 0, 1);
        else
            a = 2/(k+2);
        end
        
        % Take step
        u = u + a*(zu - u);
        v = v + a*(zv - v);
        
        %fprintf(1, '%e %e %e\n', min(u), max(0, max(u) - D), sum(u) - 1);
        %fprintf(1, '%e %e %e\n', min(v), max(0, max(v) - D), sum(v) - 1);
        assert(min(u) >= -1e-6 && max(u) <= D + 1e-6 && abs(sum(u) - 1) < 1e-6);
        assert(min(v) >= -1e-6 && max(v) <= D + 1e-6 && abs(sum(v) - 1) < 1e-6);
        
        % Calculate Lan's bound
        lan_bnd_sum_term = k / sqrt(k) * diamX + k*muk + sqrt(n)/muk*(diamX^2);
        lan_bnd_sum = lan_bnd_sum + lan_bnd_sum_term;
        lan_bnd = 2*M / (k * (k+1)) * lan_bnd_sum;
        
        lan_bnd2 = 4 * (1 + 2*n^0.25) * M * diamX /  (3 * sqrt(k));
        
        obj = norm(A*u - B*v, Inf);
        if mod(k, opts.print_interval) == 1 || opts.print_interval == 1
            [true_gu, true_gv] = gradient(A, B, u, v);
            graderr = norm([gu - true_gu; gv - true_gv]);
            
            fprintf(1, '%d %g %g %g %g %g %d %g\n', k, obj, obj - dualobj, lan_bnd2, graderr, a, nnz(grad), norm(grad));
        end
        stat.objective(k) = obj;
        stat.gap(k) = obj - dualobj;
        stat.active_indices(k) = nnz(grad);
        stat.nonzeros(k) = sum(abs([u; v]) > 1e-8);
        stat.step_nonzeros(k) = sum(abs([zu; zv]) > 1e-8);
        stat.iter_bound(k) = lan_bnd2;
  
        k = k + 1;
    end
    
    alpha = zeros(n, 1);
    alpha(y >= 0) = u;
    alpha(y < 0) = v;
    
    [w,b] =  norm1svm_primal_from_dual(X, y, alpha);
    
    % Wrap up status
    k = k - 1;
    stat.iters = k;
    stat.converged = stat.gap(end) < opts.converge_tol;
    stat.clock = toc;
end
