function H = NSDpathLinear(P,G,IC)
%computes SD paths in special case where path is linear
    a = G{2}(IC);
    hDir = 1i/a; 
    H = IC(1) + hDir*P;
end
