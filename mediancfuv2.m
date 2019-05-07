nsamples = 1e4;
d = 2;

violu = {};
violv = {};

lhs = @(u, v) 1/norm(v) * (norm(u) * norm(v) - dot(u,v));
rhs = @(u, v) norm(u-v);

for i=1:nsamples
    u = rand(d, 1) * exp(rand()*10-5);
    v = rand(d, 1) * exp(rand()*10-5);
    
    if lhs(u, v) > rhs(u, v)
        %fprintf(1, 'violation %g > %g\n', lhs, rhs);
        
        violu = [violu, u];
        violv = [violv, v];
    end
end


nviol = length(violu);
nu = zeros(1, nviol);
nv = zeros(1, nviol);
viollhs = zeros(1, nviol);
violrhs = zeros(1, nviol);

for i=1:nviol
    u = violu{i};
    v = violv{i};
    
    nu(i) = norm(u, 1);
    nv(i) = norm(v, 1);

    viollhs(i) = lhs(u, v);
    violrhs(i) = rhs(u, v);
end

%[y, i] = min(abs(nu ./ nv - 1));
if nviol ~= 0
    [y, i] = max(viollhs - violrhs);
    u = violu{i};
    v = violv{i};
end

fprintf('Violation ratio: %g\n', nviol / nsamples);
%plot([0 u(1)], [0 u(2)], [0 v(1)], [0 v(2)])