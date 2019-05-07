p = rand(5,1);
Ax = rand(5,1);
As = rand(5,1);
%As = 2*Ax - p;

d = (Ax - p) / norm(Ax - p);
Ay = @(a) Ax + a * (As - Ax);


u = @(a) Ay(a) - p;
v = Ax - p;

Cf = @(a) 1/a^2 * (norm(Ay(a) - p) - norm(Ax - p) - dot(Ay(a) - Ax, d));
%Cf2 = @(a) 1/a^2 * (norm(Ax - p + a * (As - Ax)) - norm(Ax - p) - a * dot(As - Ax, d));
%Cf2 = @(a) 1/a^2 * (2*a*norm(As - Ax) - a * dot(As - Ax, d));
%Cf2term = @(a) norm(Ay(a) - p) * (1 - (dot(Ay(a) - p, Ax - p) / (norm(Ay(a) - p) * norm(Ax - p))));
Cf2term = @(a) 1/norm(v) * (norm(u(a)) * norm(v) - dot(u(a),v));
Cf2 = @(a) 1/a^2 * Cf2term(a);

nsamples = 120;
alphasamp = zeros(nsamples, 1);
Cfsamp = zeros(nsamples, 1);
Cfsamp2 = zeros(nsamples, 1);

for i=1:nsamples
    alphasamp(i) = exp(-i/10);
    Cfsamp(i) = Cf(alphasamp(i));
    Cfsamp2(i) = Cf2(alphasamp(i));
end

semilogx(alphasamp, Cfsamp, alphasamp, Cfsamp2);
legend('Cf', 'bound')

max(abs(Cfsamp2 - Cfsamp))
max(Cfsamp2 - Cfsamp)