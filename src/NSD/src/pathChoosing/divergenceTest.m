function [divergeTF,nFinal] = divergenceTest(Z,g)
%check if SD path accidentally diverges onto a path of ascent, which can happen
%for larger #pts.
    N=length(Z);
    nFinal=N;
    divergeTF=false;
    smallErr=.1; %allow small increases in value
    for n=1:(N-1)
        if abs(exp(-imag(g(Z(n))))) + smallErr < abs(exp(-imag(g(Z(n+1)))))
            nFinal=n;
            divergeTF=true;
            warning('A steepest descent path has diverged, and will be cut short');
            return;
        end
    end
end

