p = rand(5,1);
Ax = rand(5,1);
As = rand(5,1);
%As = 2*Ax - p;

diam = norm(As - Ax);
p = p/diam;
Ax = Ax/diam;
As = As/diam;


d = (Ax - p) / norm(Ax - p);
Ay = @(a) Ax + a * (As - Ax);

u = @(a) Ay(a) - p;
v = Ax - p;

nsamples = 100;
alphasamp = zeros(nsamples, 1);
exprsamp = zeros(nsamples, 1);

for i=1:nsamples
    %alphasamp(i) = exp(-i/10);
    alphasamp(i) = (i-1)/(nsamples - 1);
    a = alphasamp(i);
    
    norm(u(a) - v) - a
    
    exprsamp(i) = 1/norm(v) * (norm(u(a)) * norm(v) - dot(u(a), v));
end

plot(alphasamp, exprsamp, alphasamp, alphasamp .^ 2)
legend('Cf', 'bound')