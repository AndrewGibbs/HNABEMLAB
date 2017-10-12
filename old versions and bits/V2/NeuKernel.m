function val=NeuKernel(kwave,R,n_dot_x_minus_source)
    %kernel of Neumann trace of single layer
    val= -(1i*kwave/4)*besselh(1,1,self.kwave*R).*n_dot_x_minus_source./R;
end