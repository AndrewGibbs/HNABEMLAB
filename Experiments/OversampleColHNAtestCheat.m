%test 'convergence' of oversampled collocation HNA BEM

clear classes;
addpath ..;
addPaths();
%ALSO NEED TO ADD SOLVER CLASS
%surpress ill conditioning warnings, as there'll be a lot of these.
warning('off','MATLAB:rankDeficientMatrix');
%set up basic problem
%define the screen
vertices=[0 0       %first vertex
          1 0];     %second vertex
%create 'edge' object for the screen
Gamma=edge(vertices);
%inident plane wave

%degree of high order standard Galerkin BEM reference solution:
standartMethodPolynomialDegree=12;

%loop parameters
P=[1:8];
K=[1024];
SampleScale=linspace(1,2,10);

kCount=0;
for kwave=K
    
    fprintf('k=%d\n',kwave);
    kCount=kCount+1;
    uinc=planeWave(kwave,[1 1]./sqrt(2));
    
    SolnName=sprintf('k%dp%d_HNA',kwave,standartMethodPolynomialDegree);
    try %see if solution was computed on a previous run
        [v,VHNAtest] = LoadSolution( SolnName );
        display('Previous solution loaded')
    catch %if not, compute it now and save it
        display('No previous solution exists');
        fprintf('Solving k=%d problem via high order HNA BEM with 2 PT/DOF \n',kwave);
        [v, VHNAtest, ~, v_Ncoeffs] = HNAColOversample( kwave,Gamma,uinc,standartMethodPolynomialDegree,2,1,1 );
        %take this solution to be the 'exact' one
        SaveSolution( VHNAtest, v_Ncoeffs, SolnName, vertices, uinc);
    end
    
    %and will use Composite Gauss for L1 error integrals:
    CG=CompGauss(kwave,30);
    [x, w]=CG.meshQuad(VHNAtest);
    L2relErr=@(v,u) (w.'*abs(v.eval(x)-u.eval(x)))/(w.'*abs(v.eval(x)));
    
    %now start the main loop over different collocations
    pCount=0;
    for p=P
        fprintf('\tp=%d\n',p);
        pCount=pCount+1;
        sCount=0;
        for s=SampleScale
            sCount=sCount+1;
            
            [v_N{pCount,kCount,sCount}, VHNA, X] = HNAColOversample( kwave,Gamma,uinc,p,s);
            DOFs(pCount,kCount,sCount)=length(VHNA.el);
%             DOFs_(pCount,kCount,sCount)=length(VHNA_.el);
%             DOFs__(pCount,kCount,sCount)=length(VHNA__.el);
            %store number of collocation points
            CPs(pCount,kCount,sCount)=length(X);
%             CPs_(pCount,kCount,sCount)=length(X_);
%             CPs__(pCount,kCount,sCount)=length(X__);
        fprintf('\t\tOversamples/element=%f:',CPs(pCount,kCount,sCount)/DOFs(pCount,kCount,sCount));
            %get the solution in each case
            %fprintf('\tNo scaling, ');
            fprintf('\tNo scaling\n');
%             fprintf('1/sqrt(k) scaling, ');
%             [v_N_{pCount,kCount,sCount}, VHNA_, X_] = HNAColOversample( kwave,Gamma,uinc,p,s,1/sqrt(kwave));
%             fprintf('1/k scaling.\n');
%             [v_N__{pCount,kCount,sCount}, VHNA__, X__] = HNAColOversample( kwave,Gamma,uinc,p,s,1/kwave);
            %now get the L1 errors, to be processed later:
            err(pCount,kCount,sCount)=L2relErr(v,v_N{pCount,kCount,sCount});
%             err_(pCount,kCount,sCount)=L2relErr(v,v_N_{pCount,kCount,sCount});
%             err__(pCount,kCount,sCount)=L2relErr(v,v_N__{pCount,kCount,sCount});
            %store DOFs
%             fprintf('\t\tErrors: %f, %f, %f\n', err(pCount,kCount,sCount), err_(pCount,kCount,sCount), err__(pCount,kCount,sCount));
            fprintf('\t\tError: %f\n', err(pCount,kCount,sCount));
            %save the output at each stage
            %save('ColTestErrorData','err','err_','err__', 'DOFs','DOFs_','DOFs__','CPs','CPs_','CPs__');
            save('ColTestErrorData','err', 'DOFs','CPs');
        end
    end
end
