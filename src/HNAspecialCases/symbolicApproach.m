syms x1 x2 yd1 yd2 P1 P2 t d1 d2 pm;

assume((yd1^4 + 2*yd1^2*yd2^2 - yd1^2 + yd2^4 - yd2^2) ==0);
%assume(- pm^2*yd1^2 - pm^2*yd2^2 + yd1^4 + 2*yd1^2*yd2^2 + yd2^4 == 0);
%assume((2*(P1*yd1^3 + P2*yd2^3 - x1*yd1^3 - x2*yd2^3 + P1*yd1*yd2^2 + P2*yd1^2*yd2 + pm^2*x1*yd1 + pm^2*x2*yd2 - x1*yd1*yd2^2 - x2*yd1^2*yd2 - P1*pm^2*yd1 - P2*pm^2*yd2)) ==0 , - pm^2*yd1^2 - pm^2*yd2^2 + yd1^4 + 2*yd1^2*yd2^2 + yd2^4 == 0);

g1 = symfun(sqrt((x1 - (P1 + t*yd1))^2+(x2 - (P2 + t*yd2))^2) + pm*t,t);
g2 = symfun(sqrt((x1 - (P1 + t*yd1))^2+(x2 - (P2 + t*yd2))^2) + (P1*d1 + P2*d2 + t*(d1*yd1 + d2*yd2)),t);
dg1 = diff(g1,t);
dg2 = diff(g2,t);
xi1 = solve(dg1 == 0,t);
xi2 = solve(dg2 == 0,t);

dg3 = symfun((yd1*(P1-x1+t*yd1) + yd2*(P2-x2+t*yd2))^2 - (P1-x1+t*yd1)^2 - (P2-x2+t*yd2)^2, t);
xi3 = solve(dg3 == 0, t);

dg4 = symfun(yd1*(P1-x1+t*yd1) + yd2*(P2-x2+t*yd2) , t);
dg5 = symfun(-pm*((P1-x1+t*yd1)^2+(P2-x2+t*yd2)^2)^.5 , t);

xi4 = solve(dg4 == dg5,t);

collect(dg3,t)

(2*(yd1^2 + yd2^2)*(yd1*(P1 - x1) + yd2*(P2 - x2)) - 2*yd2*(P2 - x2) - 2*yd1*(P1 - x1))*t + (yd1*(P1 - x1) + yd2*(P2 - x2))^2 - (P1 - x1)^2 - (P2 - x2)^2

r = symfun(sqrt((x1 - (P1 + t*yd1))^2+(x2 - (P2 + t*yd2))^2),t);
solve(r==0,t)