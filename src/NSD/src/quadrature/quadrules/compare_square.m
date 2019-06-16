% function compare_square

% Choose T larger to make the integrand simpler to approximate
T = 10;

f = @(x,y,z,t) cos( (x+y+z+t) / T);
Z = -16*sin(1/T)^2*T^4*cos(1/T)^2+16*sin(1/T)^2*T^4

% Test case: quadratic integrand
% f = @(x,y,z,t) (x+y+z+t)^2;
% Z = 64/3


[x1,w1] = encyclop_cub_square(3, 4, 2);
z = 0;
N = size(x1,1);
u = zeros(N,1);
for i=1:size(x1,1)
    f2 = @(x,y) f(x1(i,1),x1(i,2),x,y);
    u(i) = apply_cub(f2, x1, w1);
end
z1 = w1'*u;
disp('Tensor product of (3,4) rule')
abs(z1-Z)/abs(Z)

[x2,w2] = encyclop_cub_fourcube(3, 8, 2);
z2 = apply_cub4d(f, x2, w2);
disp('Four-cube rule (3,8)')
abs(z2-Z)/abs(Z)

% This is the tensor-product version of the first rule!
[x3,w3] = encyclop_cub_fourcube(3, 16);
z3 = apply_cub4d(f, x3, w3);
disp('Four-cube rule (3,16)')
abs(z3-Z)/abs(Z)


[x4,w4] = encyclop_cub_fourcube(5, 24);
z4 = apply_cub4d(f, x4, w4);
disp('Four-cube rule (5,24)')
abs(z4-Z)/abs(Z)


[x5,w5] = encyclop_cub_square(4, 6, 2);
z = 0;
N = size(x5,1);
u = zeros(N,1);
for i=1:N
    f2 = @(x,y) f(x5(i,1),x5(i,2),x,y);
    u(i) = apply_cub(f2, x5, w5);
end
z5 = w5'*u;
disp('Tensor-product of (4,6) rule')
abs(z5-Z)/abs(Z)


[x6,w6] = encyclop_cub_fourcube(7, 57);
z6 = apply_cub4d(f, x6, w6);
disp('Four-cube rule (7,57)')
abs(z6-Z)/abs(Z)

[x7,w7] = encyclop_cub_fourcube(3, 8, 1);
z7 = apply_cub4d(f, x7, w7);
disp('Four-cube rule (3,8), Stroud version')
abs(z7-Z)/abs(Z)

[x8,w8] = encyclop_cub_fourcube(3, 9);
z8 = apply_cub4d(f, x8, w8);
disp('Four-cube rule (3,9)')
abs(z7-Z)/abs(Z)
