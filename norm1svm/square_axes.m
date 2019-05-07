function square_axes(h)
    if ~exist('h', 'var')
        h = gca;
    end
    
    v = axis(h);
    sz = v([2 4]) - v([1 3]);
    
    if sz(1) < sz(2)
        grow = 0;
        from = sz(1);
        to = sz(2);
    else
        grow = 2;
        from = sz(2);
        to = sz(1);
    end
    by = to - from;
    
    v(grow + 1) = v(grow + 1) - by/2;
    v(grow + 2) = v(grow + 2) + by/2;
    axis(h, v);
end