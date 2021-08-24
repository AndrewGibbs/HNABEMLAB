z = 1i*linspace(100,2000,10000);
plot(imag(z),imag(besselh_0_1_nonosc_large_imag(z)),imag(z),imag(besselh(0,1,z)./exp(1i*z)),'-.')
legend('Combined approach','default matlab','location','northwest')
xlabel('Im(z)')
ylabel('Im[H_0^1(z)/exp(iz)]')
title('Region in complex plane where Matlab besselh stops working');