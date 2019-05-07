
a = 2;
p = 1/6;

k = 0:100;

%a = @(k) 2./(k+2);
%a = @(k) exp(-k);
a = @(k) (2./(k+2)).^(1/3);

lhs = a(k).^(1/3) - 1/2 * a(k).^(4/3);
rhs = a(k+1);

viol = find(lhs > rhs);
fprintf(1, 'Violated for:\n');
fprintf(1, 'kmin = %d\n', min(viol));
fprintf(1, 'kmax = %d\n', max(viol))