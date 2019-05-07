function [statp, statd, statfw] = run_dataset(dataset, C, D, fwopts, outpath)
    data_path = '../../data/';
    out_dir = '../../doc/nips2015/fig/svm/';
    
    % Load dataset
    dataset_path = [data_path, dataset, '.mat'];
    s = load(dataset_path);
    
    X = s.X;
    [d n] = size(X);
    fprintf(1, 'n = %d, d = %d\n', n, d);
    fprintf(1, 'Min feasible D = %g\n', 1/n);
    
    l = sort(unique(s.y));
    assert(length(l) == 2);
    y = zeros(n, 1);
    y(s.y == l(1)) = -1;
    y(s.y == l(2)) = 1;
    
    % Run solvers
    [w, b, statp] = norm1svm_primal(X, y, C);
    fprintf('Primal nonzeros = %d\n', nnz(w));
    
    [alpha, wd, bd, statd] = norm1svm_dual(X, y, D);
    fprintf(1, 'CPLEX dual objective: %g\n', statd.obj);
    fprintf(1, 'CPLEX support vectors: %d\n', sum(abs(alpha > 1e-6)));
    
    if ~fwopts.smooth
        [alphafw, wfw, bfw, statfw] = norm1svm_frank_wolfe(X, y, D, fwopts);
    else
        [alphafw, wfw, bfw, statfw] = norm1svm_sample_smoothing(X, y, D, fwopts);
    end
    
    fprintf(1, 'CPLEX dual objective: %g\n', statd.obj);
    fprintf(1, 'Frank-Wolfe dual objective: %g (%g)\n', statfw.objective(end), statfw.objective(end) / statd.obj);
    
    % Write out FW progress
    outpath = [out_dir, outpath, '-optim.dat'];
    fprintf(1, 'Writing optimization progress to "%s"\n', outpath);
    outfp = fopen(outpath, 'w+');
    fprintf(outfp, '# k\tobj\tgap\tactive indices\tnonzeros\tstep nonzeros\tbound\n');
    for k=1:statfw.iters
        fprintf(outfp, '%d\t%g\t%g\t%d\t%d\t%d\t%g\n', ...
            k, statfw.objective(k), statfw.gap(k), ...
            statfw.active_indices(k), statfw.nonzeros(k), ...
            statfw.step_nonzeros(k), statfw.iter_bound(k));
    end
    fclose(outfp);
end
