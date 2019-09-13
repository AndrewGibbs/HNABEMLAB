k = 10000;
xi = -.0001;
T = .1;
alphaPhase = 2;
Npts = 15;

%some prefactor with log singularity which localises in k
F = @(s) besselh(0,1,k*(s-xi));

%get exact value of integral
Iml = integral(@(x) F(x).*exp(1i*alphaPhase*k*x), 0, T, 'RelTol',1e-14);
singularity.position = xi;
singularity.blowUpType = 'log';
singularity.distFun = @(x) abs(x-xi);
[x, w] = NonOsc45(0,T,k,20,@(x) alphaPhase*x,singularity);
Igq = w.'*F(x);

[x_,w_] = LinearNSD(T,k,alphaPhase,xi,Npts);
I_ = w_.'*F(x_);

%print relative error:
fprintf('\ngq err: %e',abs(Igq-I_)/abs(Igq));
fprintf('\nml err: %e',abs(Iml-I_)/abs(Iml));