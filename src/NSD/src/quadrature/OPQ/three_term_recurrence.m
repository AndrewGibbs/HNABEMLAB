function z = three_term_recurrence(ab, x)

n = size(ab,1);
z = zeros(n+1,1);

z(1) = 0;
z(2) = 1;

for i = 3:n+1
    z(i) = (x-ab(i-2,1)) * z(i-1) - ab(i-2,2)*z(i-2);
end

z = z(2:end);
