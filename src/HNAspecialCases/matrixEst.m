function bound = matrixEst(L,k)
%an upper bound of the matrix entries for single layer collocation, given
%the side-length and wavenumber
    bound = sqrt(5*2.1^2/(8*log(2)*k)*log(2+k*L));
end

