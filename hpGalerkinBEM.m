


classdef hpGalerkinBEM < solver

    % ** Note: the class name above should match the filename and the name
    % of the constructor below **
    
    properties
        Gamma
        %basis bits:
        pMax
        h
        nLayers
        sigmaGrad=0.15
        V
        %quad routine:
        CG
        %discrete system bits
        GalerkinMatrix
        coeffs
        v_N
    end
    
    methods
        

        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = hpGalerkinBEM(kwave,incidentField, Gamma,p, h)

            
            % call parent constructor
            self = self@solver(kwave,incidentField);
            
            self.pMax=p;
            self.h=h;
            self.nLayers=2*(1+p);
            self.Gamma=Gamma;
            %create basis
            self.V=hpStandardBasis(self.Gamma, self.pMax, self.h, self.nLayers, self.sigmaGrad);
            %create integrator ojcect to do Galerkin integrals:
            self.CG=CompGauss(self.kwave,30);
            
        end
                 
        %===============================================================
        % these methods must be provided 
        %===============================================================

        %-----------------------------------------
        % setup
        %-----------------------------------------
        
        % This methods sets up your solver, eg assembles discretisation
        % matrices etc
        
        function setup(self)                      

            %define the single layer 'operator' object
            S=SingleLayer(self.kwave,self.Gamma);
            
            %initialise Galerkin matrix:
            GalerkinMatrix=zeros(self.V.numEls);
            for n=1:self.V.numEls
                for m=1:self.V.numEls
                    GalerkinMatrixIntegrals{n,m}=L2(S*self.V.el(m),self.V.el(n));
                    self.GalerkinMatrix(n,m)=self.CG.eval(GalerkinMatrixIntegrals{n,m});
                    %GalerkinMatrixGood(n,m)=CG_.eval(GalerkinMatrixIntegrals{n,m});
                end
                fprintf('%d / %d\n',n,self.V.numEls);
            end
            
        end
        
        %-----------------------------------------
        % solve
        %-----------------------------------------
        
        % This method solves the scattering problem for every right hand
        % side specified in the self.incidentField cell array.
        
        function solve(self)
            
            for m=1:self.numIncidentField
            %with RHS data
            f{m}=DirichletData(self.incidentField{m},self.Gamma);
                for n=1:self.V.numEls
                    GalerkinRHSIntegrals1{n,m}=L2(f{m},self.V.el(n));
                    GalerkinRHS(n,m)=self.CG.eval(GalerkinRHSIntegrals1{n});
                    fprintf('%d / %d\n',n,self.V.numEls);
                end
            end
            % ** recommend you check here that the setup method has been
            % run **
            
            self.coeffs=self.GalerkinMatrix\GalerkinRHS;     
            for m=1:self.numIncidentField
                self.v_N{m}=Projection(self.coeffs,self.V);
            end
            
        end
        
        function val = NeuTrace(self,points,index)
            if nargin <= 2
                index=1;
            end
            val=self.v_N{index}.eval(points);
        end
        
        %-----------------------------------------
        % get far field
        %-----------------------------------------
        
        % This method should compute the far field for the incident fields
        % self.incidentField{k} for each k in the array index. The return
        % value val should be an array with the column val(:,j) containing
        % the far field for self.incidentField{index(j)}.
        
        function val = getFarField(self,points,index)
            
            % ** your code goes here **            

        end
        
        %===============================================================
        % you may provide other methods required to implement your solver
        % or help the user
        %===============================================================

    end % end methods
    
end