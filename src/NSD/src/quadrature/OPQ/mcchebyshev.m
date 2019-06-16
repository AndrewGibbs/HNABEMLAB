% MCCHEBYSHEV Multiple-component discretized Chebyshev algorithm.
%
%    This is an alternative to the routine chebyshev. It performs
%    a sequence of discretizations of the given weight function
%    (or measure), each discretization being used to compute
%    discretized modified moments, which in turn are passed
%    as input to the modified Chebyshev algorithmi in order
%    to produce approximations to the desired recurrence
%    coefficients. The fineness of the discretization is
%    characterized by a discretization parameter M. The support 
%    of the continuous part of the weight function is decomposed 
%    into a given number mc of subintervals (some or all of which 
%    may be identical). The routine then applies to each 
%    subinterval an M-point quadrature rule to discretize the
%    weight function on that subinterval. The discrete part of the
%    weight function (if there is any) is added on to the
%    discretized continuous weight function. The sequence of 
%    discretizations, if chosen judiciously, leads to convergence
%    of the recurrence coefficients for the discretized measures 
%    to those of the given measure. If convergence to within a
%    prescribed accuracy eps0 occurs before M reaches its maximum
%    allowed value Mmax, then the value of M that yields 
%    convergence is output as Mcap, and so is the number of 
%    iterations, kount. If there is no convergence, the routine 
%    displays the message "Mcap exceeds Mmax in mccheb" prior to
%    exiting.
%
%    The details of the discretization are to be specified prior 
%    to calling the procedure. They are embodied in the following 
%    global parameters:
%
%    mc     = the number of component intervals
%    mp     = the number of points in the discrete part of the 
%             measure (mp=0 if there is none)
%    iq     = a parameter to be set equal to 1, if the user
%             provides his or her own quadrature routine, and 
%             different from 1 otherwise
%    idelta = a parameter whose default value is 1, but is 
%             preferably set equal to 2, if iq=1 and the user 
%             provides Gauss-type quadrature routines
% 
%    The component intervals have to be specified (in the order 
%    left to right) by a global mcx2 array AB=[[a1 b1];[a2 b2];
%    ...;[amc bmc]], where for infinite extreme intervals a1=-Inf
%    resp. bmc=Inf. The discrete spectrum (if mp>0) is similarly
%    specified by a global mpx2 array DM=[[x1 y1];[x2 y2];...;
%    [xmp ymp]] containing the abscissae and jumps.
%
%    If the user provides his or her own quadrature routine 
%    "quadown", the routine mcdis must be called with the input 
%    parameter "quad" replaced by "@quadown", otherwise with 
%    "quad" replaced by "@quadgp", a general-purpose routine
%    provided in the package. The quadrature routine must have 
%    the form
%
%                             function xw=quad(M,i)
%
%    where M is the number of nodes and i identifies the interval
%    to which the routine is to be applied.
%
%    The routine mcchebyshev also applies to measures given
%    originally in multi-component form.
%
function [ab,Mcap,kount]=mcchebyshev(N,eps0,quad,Mmax)
global mc mp iq idelta DM uv AB 
global om2 abm
if N<1, error('Input variable N out of range'), end
iMcap=1; kount=-1;
ab(:,2)=zeros(N,1); b=ones(N,1);
Mcap=floor((2*N-1)/idelta);
while any(abs(ab(:,2)-b)>eps0*abs(ab(:,2)))
  b=ab(:,2);
  kount=kount+1;
  if kount>1, iMcap=2^(floor(kount/5))*N; end
  Mcap=Mcap+iMcap;
  if Mcap>Mmax, error('Mcap exceeds Mmax in mccheb'), end
  mtMcap=mc*Mcap;
  if iq~=1, uv=fejer(Mcap); end
  for i=1:mc
    im1tn=(i-1)*Mcap;
    xw=feval(quad,Mcap,i);
    xwm(im1tn+1:im1tn+Mcap,:)=xw;
  end
  if mp~=0, xwm(mtMcap+1:mtMcap+mp,:)=DM; end
  for k=1:2*N
    p1=zeros(mtMcap+mp,1);
    p=ones(mtMcap+mp,1);
    if k>1
      for l=1:k-1
        pm1=p1;
        p1=p;
        p=(xwm(:,1)-abm(l,1)).*p1-abm(l,2)*pm1;
      end
    end
    mom(k)=xwm(:,2)'*p;
  end
  ab=chebyshev(N,mom,abm);
end
