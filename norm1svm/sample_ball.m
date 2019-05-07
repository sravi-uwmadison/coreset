function eta = sample_ball(n, r)
    % eta = sample_ball(n)
    % eta = sample_ball(n, r)
    %
    %    Sample uniformly from an n-dimensional ball around the
    %    original of radius r.  Default value of r is 1.
    %

    if ~exist('r', 'var')
        r = 1;
    end
    
    % Rejection sampling approach: will not scale to large # of dimensions
    %eta = rand(n, 1);
    %while norm(eta) > 1:
    %    eta = rand(n, 1);
    %    continue
    %end
    
    % Sample on unit sphere by scaling an isotropic multivariate
    % normal
    eta = normrnd(0, 1, n, 1);
    eta = eta / norm(eta);
    
    % Sample the "radius" term so that we can also cover the
    % interior of the unit ball
    R = betarnd(n, 1);
    eta = eta * R;
    
    assert(norm(eta) <= 1 + 1e-8)
    eta = eta * r;
end