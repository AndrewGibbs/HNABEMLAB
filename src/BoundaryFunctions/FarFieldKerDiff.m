function [ker_vals, ker_vals_all_derivs] = FarFieldKerDiff( theta, y, k, n)
%Function to compute the nth derivative of the far field kernel. Symbolic
%matab takes so long to do this for a general n, but I think it can be done
%quickly by switching to a Fourier basis.
   %y and theta are the inputs, in vertical vector form (y is a 2*N matrix), whilst dn is the
   %number of derivatives to apply to the far field kernel. Wavenumber is
   %k, as usual.
   
   %the original far field kernel is:
   ker=@(theta,y) exp(-1i*k*(y(:,1)*cos(theta.')+y(:,2)*sin(theta.')));
   ker_vals_all_derivs=cell(n+1,1);
   if nargin<4
       n=0;
   end
   
   %if n==0
       %simple case, no derivatives of far field kernel
       %ker_vals=ker(theta,y);
       %Compute n=0 derivative
       ker_vals_all_derivs{1}=ker(theta,y);
   %else
    if n>0
       %STEP 1: switch to exponential Fourier basis
        [h_y,h_w_test]=size(y);
        if h_w_test~=2
            error('y needs to be two column vectors');
        end
        c=[y(:,2)-1i*y(:,1) zeros(h_y,1) y(:,2)+1i*y(:,1)]*-1i*k/2;

       %Compute n=1 derivative
        Fprod=c;
        e_ells_thetas=exp(1i*theta*(-1:1)).';
        ker_vals_all_derivs{2} = ker(theta,y) .* (Fprod * e_ells_thetas) ;
    end
        
       %Compute n>1 derivatives
        if n>1
            for n_=2:n
              %the effect on the coefficients by a differentiation in theta
               DerF=1i*[zeros(h_y,1) repmat((-(n_-1):(n_-1)),h_y,1) .*Fprod    zeros(h_y,1)];
              %the effect on the coefficients by a multiplication with f_1
               fF=repmat(c(:,1),1,2*n_+1).*[Fprod zeros(h_y,2)] + repmat(c(:,3),1,2*n_+1).*[zeros(h_y,2) Fprod];
               Fprod=DerF+fF;

                %STEP 3: Create matrix of exponential basis terms exp(i\ell\theta)
                e_ells_thetas=exp(1i*theta*(-n_:n_)).';
                
                %STEP 4: Mash it all together
                ker_vals_all_derivs{n_+1} = ker(theta,y) .* (Fprod * e_ells_thetas) ;
            end
        end

   %end
    ker_vals=ker_vals_all_derivs{n+1};
end