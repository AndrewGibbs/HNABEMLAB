% *** make sure ChebFun and PathFinder are added to path before running this
% code***
clear classes;
a = 0.1;
b = 1;
freq = 1000;
N = 20;

%mimic the type of function we're likely to get in the HNA code:
SING = a-0.01;
amp = @(x) log(freq*(x-SING));
fSing = singularity(SING, 'log', @(r) abs(r-SING));
y = -SING; theta = pi/2;
%phase = @(x) sqrt(x.^2 + y.^2 - 2*x.*y.*cos(theta));
phase = @(x) (x-2).^2;

%now ChebFun it:
AMP = chebfun(amp,[a b]);
chebPhase = chebfun(phase,[a b]);
DchebPhase = diff(chebPhase);
DDchebPhase = diff(diff((chebPhase)));
G = {@(x) chebPhase(x), @(x) DchebPhase(x)};    % @(x) DDchebPhase(x)};

%now compute SD path, and evaluate integrals:
[ z, w ] = PathFinder( a, b, freq, N, G, 'stationary points', [], 'order', [],'fsingularities',fSing,'minoscs',20);

I = integral(@(x) amp(x).*exp(1i*freq*phase(x)),a,b,'RelTol',1e-16);
I_N = I;%(F(z).*exp(1i*freq*chebPhase(z)));
I_Ngrad = PathFinderChebWrapGrad(SING, a, b, freq, N, amp, G);
err1 = abs(I-I_N) /I
err2 = abs(I-I_Ngrad) /I
plot(z,'x');