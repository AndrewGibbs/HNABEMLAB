function [theta, FFcoarse] = getFFtroughs(theta, FFcoarse, logAbsFFhandle)
    tol = 1e-12;
    %approximate the sign of the derivatives:
    FFcoarse = [FFcoarse(end); FFcoarse(:); FFcoarse(1);];
    theta = [theta(end); theta(:); theta(1);];
    Dsgn = sign(FFcoarse(2:end) - FFcoarse(1:(end-1)));
    
    %isolate the dodgy indices
    oneToN = 1:length(Dsgn);
    desiredSignChanges = oneToN((Dsgn(1:(end-1))==-1)&(Dsgn(2:end)==1))+1;
    
    %adjust them to the exact trough point
    for n=desiredSignChanges
        if theta(n-1) < theta(n+1)
            [theta(n), FFcoarse(n)] = ...
            fminbnd(logAbsFFhandle,theta(n-1),theta(n+1),optimset('TolX',tol,'Display','off'));
        else
            [theta(n), FFcoarse(n)] = ...
            fminbnd(logAbsFFhandle,theta(n-1)-2*pi,theta(n+1),optimset('TolX',tol,'Display','off'));
            if theta(n)<0
                theta(n) = theta(n) + 2*pi;
            end
        end
    end
    theta = theta(2:(end-1));
    FFcoarse = FFcoarse(2:(end-1));
end