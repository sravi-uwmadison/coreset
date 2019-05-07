
f = @(x) abs(x);
y = @(x,s,a) x + a*(s - x);

cfterm = @(x, s, a, d) (f(y(x,s,a)) - f(x) - d'*(y(x,s,a)-x))/a^2;

eps = 1e-4;
x = -eps/2;
d = -1;
s = 1 - eps/2;
a = eps;

y(x,s,a)
f(y(x,s,a))
f(x)
d'*(y(x,s,a)-x)
f(y(x,s,a)) - f(x) - d'*(y(x,s,a)-x)
cfterm(x, s, a, d)