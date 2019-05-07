rng(0);

n = 100;
d = 2;

mu1 = [rand(1, d)/5];
mu2 = [rand(1, d)/5];

active = ceil(rand()*d/10)
active = randsample(d, active, 0);
mu1(active) = 0;
mu2(active) = 4;

mu1 = [0, 0.5];
mu2 = [4, -0.5];

sigma = 1e-1;

%C = 0.2;
C = 1000;

% Generate labels
y = 2*(rand(n, 1) < 0.5) -1;
nneg = sum(y == -1);
npos = sum(y == 1);
assert(n == nneg + npos);

% Generate points
X = zeros(d, n);
X(:, y == -1) = mvnrnd(mu1, eye(d) * sigma, nneg)';
X(:, y == 1) = mvnrnd(mu2, eye(d) * sigma, npos)';

% Solve svm
[w, b, statp] = norm1svm_primal(X, y, C);
%figure; svmplot(X, y, w, b);

%D = 1/10;
D = 1;
[alpha, wd, bd, statd] = norm1svm_dual(X, y, D);
fprintf(1, 'True min: %g\n', norm(X * (y .* alpha), Inf))

opts.max_iters = 100;
[alphafw, wfw, bfw, statfw] = norm1svm_frank_wolfe(X, y, D, opts);
%fprintf(1, 'F-W min: %g\n', norm(X * (y .* alphafw), Inf))

opts.print_interval = 1;
[alphasm, wsm, bsm, statsm] = norm1svm_sample_smoothing(X, y, D, opts);
%fprintf(1, 'F-W min: %g\n', norm(X * (y .* alphafw), Inf))

fprintf(1, 'CPLEX:\n');
fprintf(1, '  Primal nonzeros: %d\n', nnz(w));
fprintf(1, '  %d support vectors\n', nnz(alpha))

fprintf(1, 'Frank-Wolfe solution:\n');
fprintf(1, '  Primal nonzeros: %d\n', nnz(wfw));
fprintf(1, '  %d support vectors\n', sum(abs(alphafw) > 1e-10));

fprintf(1, 'Frank-Wolfe w/ Randomized Smoothing solution:\n');
fprintf(1, '  Primal nonzeros: %d\n', nnz(wsm));
fprintf(1, '  %d support vectors\n', sum(abs(alphasm) > 1e-10));

svmplot(X, y, wfw, bfw);
