fwopts.max_iters = 1000;
fwopts.line_search = 0;
fwopts.print_interval = 1;
fwopts.converge_tol = -Inf;

fwopts.smooth = 1;
suffix = '-smooth';

% Seperable?
fwopts.max_iters = 1000;
[statp, statd, statfw] = run_dataset('leukemia/leu', Inf, 1, fwopts, ['leukemia', suffix]);

fwopts.max_iters = 1000;
[statp, statd, statfw] = run_dataset('colon-cancer/colon-cancer', Inf, 1, fwopts, ['colon-cancer', suffix]);

% Good?
fwopts.max_iters = 1000;
[statp, statd, statfw] = run_dataset('ionosphere/ionosphere_scale', 10, 1e-2, fwopts, ['ionosphere', suffix]);

%fwopts.max_iters = 1000;
%[statp, statd, statfw] = run_dataset('w8a/w8a', 10, 8e-4, fwopts, ['w8a', suffix]);
%[statp, statd, statfw] = run_dataset('breast-cancer/breast-cancer', 10, 1e-2, fwopts, ['breast-cancer', suffix]);

% Unknown
%[statp, statd, statfw] = run_dataset('real-sim/real-sim', 10, 8e-4);
%[statp, statd, statfw] = run_dataset('news20.binary/news20.binary', Inf, 1);

% Ugly
%[statp, statd, statfw] = run_dataset('diabetes/diabetes_scale', 10, 5e-3);

% Bad
%[statp, statd, statfw] = run_dataset('a9a/a9a', 10, 2e-4);
%[statp, statd, statfw] = run_dataset('mushrooms/mushrooms', 10, 1e-3);
%[statp, statd, statfw] = run_dataset('covtype/covtype.libsvm.binary', 10, 1e-3);
