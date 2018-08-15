function I = LevinQuick(a,b,f,G,k)
%Quick Levin method, for use in case without stationary points
    %define the Levin ODE:
    Ly = @(x,y) f(x) -  1i*k*G{2}(x)*y;
    %solve it:
    [~,y] = ode45(Ly,[a b],0,odeset('RelTol',1E-8));
    %solution is anti-derivative, so use this to estiamate integral:
    I = y(end)*exp(1i*k*G{1}(b));%-y(1)*exp(1i*k*G{1}(a));
end