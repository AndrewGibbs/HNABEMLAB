function bound = rhsEst(L,k)
%an upper bound of the rhs entries for single layer collocation, given
%the side-length and wavenumber
    bound = L*2*k*sqrt(5*2.1^2/(8*log(2)*k)*log(2+k*L))  +  1;
end