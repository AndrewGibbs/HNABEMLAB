function [Zout, dhout]= pathAdjust(G, P, Z, h0, order, thresh)
%use householder methods, as close to a stationary point

    F = {@(w, q) G{1}(w)-(1i*q^order +G{1}(h0)), @(w) G{2}(w), @(w) G{3}(w)};
    Zout = NaN(size(Z));
    dhout = Zout;
    maxNumItrs = 100;
    
    for j = 1:length(Z)
        z = Z(j);
        p = P(j);
        initErr = abs(F{1}(z,p));
        err = initErr;
        itrCount = 0;
        while thresh < err
            %if order == 1
                z = z - F{1}(z,p)/F{2}(z);
                %use Halley's method instead:
                
                 %z = z - 2*F{1}(z,p).*F{2}(z)./(2*(F{2}(z)).^2 - F{1}(z,p)*F{3}(z) );
    %         else
    %             error('Need to code a higher order Householder method here');
    %         end
            err = abs(F{1}(z,p));
            itrCount = itrCount + 1;
            if  itrCount >= maxNumItrs
                if err>=initErr
                    Zout(j) = Z(j);
                end
                break;
            end
        end
        dhout(j) = order*1i*p^(order-1)/G{2}(z);
        Zout(j) = z;
        clear z;
    end
end