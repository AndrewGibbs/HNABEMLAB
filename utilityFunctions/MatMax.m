function [ max_el, indices] = MatMax( A )
%finds maximum entry of a matrix because matlab can't fucking do it
    sizeA=size(A);
    M=sizeA(1); N=sizeA(2);
    max_el=-inf;
    for m=1:M
        for n=1:N
            if A(m,n)>max_el
                max_el=A(m,n);
                indices=[m,n];
            end
        end    
    end
end

