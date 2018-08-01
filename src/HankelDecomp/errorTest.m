% theta = linspace(0,2*pi);
% r = linspace(0,1000);
% z = r.'*exp(1i*theta);

x = linspace(-5000,5000);
z = x + 1i*x.';

Hz = besselh(0,1,z);
HzDC = besselhDecomp(0,1,z);

relErr = abs(Hz - HzDC)./abs(Hz);

figure(1); imagesc(x,x,log10(relErr));

figure(2); imagesc(isinf(Hz));