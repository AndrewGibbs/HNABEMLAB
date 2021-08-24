function a = Hankel_asym_coeffs(k,v)
    if k==0
        a = 1;
    elseif v==0 % some presets, for commonly used stuff
        switch k
            case 1
                a = -0.125;
            case 2
                a = 0.0703125;
            case 3
                a = -0.0732421875;
            case 4
                a = 0.112152099609375;
            case 5
                a = -0.227108001708984;
            case 6
                a = 0.572501420974731;
            case 7
                a = -1.727727502584457;
            case 8
                a = 6.074042001273483;
            case 9
                a = -24.380529699556064;
            case 10
                a = 1.100171402692467e+02;
            case 11
                a = -5.513358961220206e+02;
            otherwise
                a = compute_a(k,v);
        end
    else
        a = compute_a(k,v);
    end
    function a_fresh = compute_a(k,v)
        a_fresh = prod(4*v^2 - (1:2:(2*k-1)).^2)/(factorial(k)*8^k);
    end
end