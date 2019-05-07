n = 10000;
d = 5;
nCfsamp = 2500;

%p = [-1, 1, 0.1 * rand(1, n)];
A = rand(d, n);

colnorms = @(A) sqrt(sum(A .^ 2, 1));
Adists = @(x) colnorms(repmat(x, 1, n) - A);

fhat = @(x) sum(Adists(x));
dfhat = @(x) sum((repmat(x, 1, n) - A) ./ repmat(Adists(x), d, 1), 2);
Cfhat = @(x,s,alpha) (fhat(x + alpha*(s-x)) - fhat(x) - alpha * (s - x)' * dfhat(x)) / alpha^2;

maxCfhat = -Inf;
for i=1:nCfsamp
    x = rand(d, 1);
    s = rand(d, 1);
    alpha = rand(1);
    
    maxCfhat = max(maxCfhat, Cfhat(x, s, alpha));
end

maxCfhat

f = @(x) fhat(A*x)/n;
df = @(x) A' * dfhat(A*x)/n;
Cf = @(x,s,alpha) (f(x + alpha*(s-x)) - f(x) - alpha * (s - x)' * df(x)) / alpha^2;

maxCf = -Inf;
for i=1:nCfsamp
    x = -log(rand(n,1));
    x = x / sum(x);
    %x = zeros(n,1);
    %x(randsample(n,1)) = 1;
    
    if rand() < 0.5
        s = -log(rand(n,1));
        s = s / sum(s);
    else
        s = zeros(n,1);
        s(randsample(n,1)) = 1;
    end
    
    if rand() < 0.3
        alpha = rand() * 1e-4;
    elseif rand() < 0.6
        alpha = 1 - rand() * 1e-4;
    else
        alpha = rand();
    end
    
    %Cf(x, s, alpha)
    maxCf = max(maxCf, Cf(x, s, alpha));
end

maxCf