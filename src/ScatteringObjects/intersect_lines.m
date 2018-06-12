function S = intersect_lines( X,X_0,Y,Y_0 )
    %Finds the values s and t such that
        % X_0 + sX = Y_0 + tY
    S=[X,Y]\(Y_0-X_0);
    %s(1)=S(1);
    S(2)=-S(2);
end