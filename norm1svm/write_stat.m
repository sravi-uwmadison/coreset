function write_stat(statfw, outpath)
    out_dir = '../../doc/nips2015/fig/svm/';

    % Write out FW progress
    outpath = [out_dir, outpath, '-optim.dat'];
    fprintf(1, 'Writing optimization progress to "%s"\n', outpath);
    outfp = fopen(outpath, 'w+');
    fprintf(outfp, '# k\tobj\tgap\tactive indices\tnonzeros\tstep nonzeros\n');
    for k=1:statfw.iters
        fprintf(outfp, '%d\t%g\t%g\t%d\t%d\t%d\n', k, statfw.objective(k), statfw.gap(k), statfw.active_indices(k), statfw.nonzeros(k), statfw.step_nonzeros(k));
    end
    fclose(outfp);
end