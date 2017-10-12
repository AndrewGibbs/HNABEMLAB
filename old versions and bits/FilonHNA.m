classdef FilonHNA < Integrator
    %adaptation of Steve's old integration routines for the screen
    
    properties
        kwave
        meshEls %mesh elements of HNAbasis
        Npts %this is Nq1...Nq5 in Steve's code
        nlx %number of layers on screen, n_l^+
        nly %number of layers on screen, n_l^-
        gptsx %vector of gridpts supporting e^{iks}
        gptsy %vector of gridpts supporting e^{-iks}
        gdiffx %interval length supporting e^{iks}
        gdiffy %interval length supporting e^{-iks}
        plx %poly degree on each element, p^+
        ply %poly degree on each element, p^-
        Nx %total number of e^{iks} basis functions
        Ny %total number of e^{-iks} basis functions
        N %total number of degrees of freedom
    end
    
    methods
        
        %constructor
        function self=FilonHNA(HNAbasis,Npts)
            if ~isa(HNAbasis,'HNAoverlappingMesh')
               error('Steves Filon routine set up for overlapping mesh HNA basis only');
            end
            
            %MAKE SURE VECTORS ARE SAME SHAPE AS STEVE'S CODE
            self.nlx=length(HNAbasis.mesh{1}.el);
            self.nly=length(HNAbasis.mesh{2}.el);
            self.gptsx=HNAbasis.mesh{1}.points;
            self.gptsy=HNAbasis.mesh{2}.points;
            self.Nx=length(HNAbasis.mesh{1}.el);
            self.Ny=length(HNAbasis.mesh{2}.el);
            for m=1:self.Nx
                self.gdiffx(m)=HNAbasis.mesh{1}.el(m).width;
            end
            for m=1:self.Ny
                self.gdiffy(m)=HNAbasis.mesh{2}.el(m).width;
            end
            
            
            self.Nx=0;
            self.Ny=0;
            self.plx=zeros(self.Nx,1);
            self.ply=zeros(self.Ny,1);
            
            self.Npts=Npts;
            self.kwave=HNAbasis.el(1).kwave;
            %loop over all mesh elements and get relevant data
            if iscell(HNAbasis.mesh)
                self.meshEls=[];
               for LR=1:length(HNAbasis.mesh) 
                    self.meshEls=[self.meshEls HNAbasis.mesh{LR}.el];
               end
            else
                self.meshEls=HNAbasis.mesh.el;
           end
            
 
           for meshEl= self.meshEls
                %for each mesh element  
                 %'self.r'...forgot what I was writing
                      
            end
        end
        
        function val=eval_LHS(integrals)
            %should be do-able just like in Steve's code, maybe a tiny
            %extra case for non-oscillatory basis elements
        end
        
        function val=eval_fRHS(integrals)
            %need a check that this is an incident wave
            
        end
        
        function val=eval_KRHS(integrals)
            %need a check that this is only a double integral, and incident
            %wave GOA
        end
    end
    
end

