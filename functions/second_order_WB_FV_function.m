%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CODE FOR THE PAPER "WELL-BALANCED FINITE VOLUME SCHEMES FOR HYDRODYNAMIC
% EQUATIONS WITH GENERAL FREE ENERGY"
% 
% AUTHOR OF THE CODE: SERGIO P. PEREZ
%
% COAUTHORS: JOSÉ A. CARRILLO, SERAFIM KALLIADASIS, CHI-WANG SHU
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FILE: second_order_WB_FV_function.m
%
% DESCRIPTION: compute the temporal derivative of the variable vector U
% with a second-oder well-balanced scheme, for a variety of free-energy
% funcitonal choices.
%
% INPUTS:
%     x-> 1D mesh with n nodes located at the centre of the finite volume 
%         cells.
%     deltax-> vector with the sizes of each of the finite volume cells
%     U-> vector of variables (density and momentum) with length 2*n.
%     pd-> power of the density in the pressure function
%         P(\rho)=nu*\rho^(pd)
%     nu-> coefficient of the density in the pressure function
%         P(\rho)=nu*\rho^(pd)
%     cik-> logical value representing the choice of including the
%         interaction kernel W(x) convoluted with the density
%     alpha-> power and coefficient of the position vector x in the kernel 
%         function  W(x)=(abs(x)^alpha)/alpha
%     cep-> logical value representing the choice of including the
%         external potential  V(x)
%     a-> coefficient of x^4 in the external potential 
%         V(x)=a*x^4+b*x^2
%     b-> coefficient of x^2 in the external potential 
%         V(x)=a*x^4+b*x^2
%     cefrHR-> logical value representing the choice of including the
%         excessive free energy of hard rods
%     M-> matrix computed to evaluate the excessive free energy of hard
%         rods, and computed only once in the main file for efficiency
%     gamma-> coefficient of the linear damping
%     cCS-> logical value representing the choice of including the
%         Cucker-Smale damping terms from collective behaviour
%     bc-> value indicating the chosen boundary conditions. If bc==1, the
%         boundary conditions are no flux, and if bc==2 they are periodic
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  



function dUdt = second_order_WB_FV_function(x,deltax,U,pd,nu,cik,alpha,cep,a,b,cefrHR,M,gamma,cCS,bc)

n=length(U)/2;

% bc==1 No flux conditions
if bc==2 % Periodic boundary conditions in both density and momentum
    U=[U(n-2);U(1:n-2);U(1);U(2*n-4);U(n-1:2*n-4);U(n-1)];
    x=[x(end);x;x(1)];
    deltax=[deltax(end);deltax;deltax(1)];
elseif bc==3 % Periodic boundary conditions in density and reflective in momentum
    U=[U(n-2);U(1:n-2);U(1);U(n-1)/U(1)*U(n-2);U(n-1:2*n-4);U(2*n-4)/U(n-2)*U(1)];
    x=[x(end);x;x(1)];
    deltax=[deltax(end);deltax;deltax(1)];
elseif bc==4 % Reflective boundary conditions in density and momentum
    U=[U(1);U(1:n-2);U(n-2);U(n-1);U(n-1:2*n-4);U(2*n-4)];
    x=[x(1);x;x(end)];
    deltax=[deltax(1);deltax;deltax(end)];
end
%--------------------------------------------------------------------------
% Evaluation of H(x,\rho): the part of the variation of the free energy 
% with respect to the density with the potentials and excessive free energy
%--------------------------------------------------------------------------


% a) Interacting potential W(x) convoluted with density

if cik==0
    Wconvrho=zeros(length(x),1);
elseif cik==1
    Wconvrho=Wconvrhofunction(x,U(1:n),alpha,deltax);
elseif cik==2
    Wconvrho=Wgaussianconvrhofunction(x,U(1:n),deltax);
end

% b) External potential V(x)

if cep==0
    V=zeros(length(x),1);
elseif cep==1
    V=a*abs(x).^4+b*abs(x).^2;
    % Other possible choices of the external potential:
    % V=-log(1+exp(-16*x.^2));
    % V=1/6*exp(-x.^2/6);
end

% c) Hard rods excessive free energy

if cefrHR==0
    HR=zeros(length(x),1);
elseif cefrHR==1
    HR=-log(1-M*U(1:n))+rot90(M,2)*(U(1:n)./(1-M*U(1:n)));
end

% d) Add all the terms to form H(x,\rho)

H=Wconvrho+V+HR;

%--------------------------------------------------------------------------
% Construction of left and right values of the variable vector U at the 
% boundaries of the cells, by means of MUSCL procedure
%--------------------------------------------------------------------------

% Reconstruction of the density
[rhoir, rhoil]=MUSCLreconstruction2(U(1:n));
   
% Reconstruction of the velocity
[uir, uil]=MUSCLreconstruction2vel(U(n+1:2*n)./U(1:n),rhoir, rhoil,U(1:n));
%[uir, uil]=MUSCLreconstruction(U(n+1:2*n)./U(1:n));
  
% Reconstruction of the source: like in the book of Bouchut
if pd>1  
%[Hir, Hil]=MUSCLreconstruction(H);  

[PiHir, PiHil]=MUSCLreconstruction2(2*U(1:n)+H);            
Hir=PiHir-pd*(rhoir).^(pd-1);
Hil=PiHil-pd*(rhoil).^(pd-1);  
      
elseif pd==1      
%[Hir, Hil]=MUSCLreconstruction(H);

[PiHir, PiHil]=MUSCLreconstruction2(log(U(1:n))+H);
[piir, piil]=MUSCLreconstruction2(log(U(1:n)));
 
Hir=PiHir-piir;
Hil=PiHil-piil;      
end

  
%--------------------------------------------------------------------------
% Construction of plus and minus values of the variable vector U at the 
% boundaries of the cells, by means of the well-balanced procedure 
%--------------------------------------------------------------------------

% Reconstruction of + and - values
  
Uiplushalfminus=zeros(2*n,1);
Uiplushalfplus=zeros(2*n,1);
  
if pd>1
      Uiplushalfplus(1:n)=((pd-1)/pd*max(((pd/(pd-1)*circshift(rhoil,-1).^(pd-1)-max(Hir,circshift(Hil,-1))./nu+circshift(Hil,-1)./nu)),0)).^(1/(pd-1));
      Uiplushalfminus(1:n)=((pd-1)/pd*max(((pd/(pd-1)*rhoir.^(pd-1)-max(Hir,circshift(Hil,-1))./nu+Hir./nu)),0)).^(1/(pd-1));
elseif pd==1
      Uiplushalfplus(1:n)=circshift(rhoil,-1).*exp((-max(Hir,circshift(Hil,-1))+circshift(Hil,-1))./nu);
      Uiplushalfminus(1:n)=rhoir.*exp((-max(Hir,circshift(Hil,-1))+Hir)./nu);
end

if bc==1 % Implement no flux conditions with bc==1
        % If bc==2 the boundary conditions are periodic
Uiplushalfplus(n)=0;
Uiplushalfminus(n)=0;
end
 
Uiplushalfplus(n+1:2*n)=Uiplushalfplus(1:n).*circshift(U(n+1:2*n)./U(1:n),-1);
Uiplushalfminus(n+1:2*n)=Uiplushalfminus(1:n).*U(n+1:2*n)./U(1:n);

%--------------------------------------------------------------------------
% Consturction of numerical flux whose inputs are the plus and minus values
% With pd>1 vacuum regions are formed and the numerical flux is kinetic
% Reference: "Kinetic formulation of conservation laws", Benoit Perthame, page 171
% With pd==1 the numerical flux is the popular local Lax-Friedrich
%--------------------------------------------------------------------------
  
  if pd>1      
      theta=(pd-1)/2;
      Aplusrho=1/sqrt(48*nu)*Uiplushalfminus(1:n).^(1-theta).*(max(zeros(n,1),Uiplushalfminus(n+1:2*n)./Uiplushalfminus(1:n)+sqrt(3*nu)*Uiplushalfminus(1:n).^theta).^2-max(zeros(n,1),Uiplushalfminus(n+1:2*n)./Uiplushalfminus(1:n)-sqrt(3*nu)*Uiplushalfminus(1:n).^theta).^2);
      Aminusrho=1/sqrt(48*nu)*Uiplushalfplus(1:n).^(1-theta).*(min(zeros(n,1),Uiplushalfplus(n+1:2*n)./Uiplushalfplus(1:n)+sqrt(3*nu)*Uiplushalfplus(1:n).^theta).^2-min(zeros(n,1),Uiplushalfplus(n+1:2*n)./Uiplushalfplus(1:n)-sqrt(3*nu)*Uiplushalfplus(1:n).^theta).^2);
      
      Aplusrhou=1/sqrt(108*nu)*Uiplushalfminus(1:n).^(1-theta).*(max(zeros(n,1),Uiplushalfminus(n+1:2*n)./Uiplushalfminus(1:n)+sqrt(3*nu)*Uiplushalfminus(1:n).^theta).^3-max(zeros(n,1),Uiplushalfminus(n+1:2*n)./Uiplushalfminus(1:n)-sqrt(3*nu)*Uiplushalfminus(1:n).^theta).^3);
      Aminusrhou=1/sqrt(108*nu)*Uiplushalfplus(1:n).^(1-theta).*(min(zeros(n,1),Uiplushalfplus(n+1:2*n)./Uiplushalfplus(1:n)+sqrt(3*nu)*Uiplushalfplus(1:n).^theta).^3-min(zeros(n,1),Uiplushalfplus(n+1:2*n)./Uiplushalfplus(1:n)-sqrt(3*nu)*Uiplushalfplus(1:n).^theta).^3);
      
      F_iplushalf=zeros(2*n,1);
      F_iplushalf(1:n)=Aplusrho+Aminusrho;
      F_iplushalf(n+1:2*n)=Aplusrhou+Aminusrhou;  
      
  elseif pd==1
      lambdamax=zeros(n,1);
      lambdamax(1:n)=max([abs(Uiplushalfplus(n+1:2*n)./Uiplushalfplus(1:n)+sqrt(nu)) abs(Uiplushalfplus(n+1:2*n)./Uiplushalfplus(1:n)-sqrt(nu)) abs(Uiplushalfminus(n+1:2*n)./Uiplushalfminus(1:n)+sqrt(nu)) abs(Uiplushalfminus(n+1:2*n)./Uiplushalfminus(1:n)-sqrt(nu))],[],2);
      
      F_iplushalf=zeros(2*n,1);
      F_iplushalf(1:n)=0.5*(Uiplushalfplus(n+1:2*n)+Uiplushalfminus(n+1:2*n))-0.5*lambdamax(1:n).*(Uiplushalfplus(1:n)-Uiplushalfminus(1:n));
      F_iplushalf(n+1:2*n)=0.5*(Uiplushalfplus(n+1:2*n).^2./Uiplushalfplus(1:n)+nu*Uiplushalfplus(1:n)+Uiplushalfminus(n+1:2*n).^2./Uiplushalfminus(1:n)+nu*Uiplushalfminus(1:n))-0.5*lambdamax(1:n).*(Uiplushalfplus(n+1:2*n)-Uiplushalfminus(n+1:2*n));      
  end
  
  F_iplushalf(isnan(F_iplushalf))=0; % remove NaN caused by no flux conditions
  F_iminushalf=zeros(2*n,1);
  F_iminushalf(1:n)=circshift(F_iplushalf(1:n),1);
  F_iminushalf(n+1:2*n)=circshift(F_iplushalf(n+1:2*n),1);
  
%--------------------------------------------------------------------------
% Construction of source term
%--------------------------------------------------------------------------  


  Sic=zeros(n,1);
  

      if pd>1
          
          Siplushalfminus=1./deltax.*(nu*Uiplushalfminus(1:n).^pd-nu*rhoir.^pd);
          
          Siminushalfplus=1./deltax.*(nu*rhoil.^pd-nu*circshift(Uiplushalfplus(1:n),+1).^pd);
          
%                      rhoilaver=max(rhoil-0.25*(Hir+Hil)+0.5*Hil,0);
%                    rhoiraver=max(rhoir-0.25*(Hir+Hil)+0.5*Hir,0);
                    
           rhoilaver=rhoil-(Hir+Hil)/2/pd+Hil/pd;
          rhoiraver=rhoir-(Hir+Hil)/2/pd+Hir/pd;
          
           Sic=1/deltax*(rhoir.^pd-rhoiraver.^pd-rhoil.^pd+rhoilaver.^pd);
%           Sic=-1/deltax*(rhoil+rhoir)/2.*(Hir-Hil);
      elseif choicepressure==1
          
          Siplushalfminus=1/deltax*(Uiplushalfminus(1:n)-rhoir);
          
          Siminushalfplus=1/deltax*(rhoil-circshift(Uiplushalfplus(1:n),+1));
          
          rhoilaver=rhoil.*exp(-0.5*(Hir+Hil)+Hil);
          rhoiraver=rhoir.*exp(-0.5*(Hir+Hil)+Hir);
          
%            rhoilaver=exp(-0.5*(Hir+Hil)+Hil+piil);
%           rhoiraver=exp(-0.5*(Hir+Hil)+Hir+piir);
%           
          Sic=1/deltax*(rhoir-rhoiraver-rhoil+rhoilaver);

          
          
      end
      
      S(n+1:2*n)=Siplushalfminus+Siminushalfplus+Sic-gamma*U(1+n:2*n);

if cCS==1
    S(n+1:2*n)=S(n+1:2*n)-CSconvfunction(x,U(1:n),U(n+1:2*n),deltax);
end

%--------------------------------------------------------------------------
% Finally, computation of the temporal derivative of the variables U
%--------------------------------------------------------------------------  

dUdt=-1./[deltax;deltax].*(F_iplushalf-F_iminushalf)+S;
  
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Auxiliary functions to compute the convolutions
%--------------------------------------------------------------------------

  
      function [ir, il]=MUSCLreconstruction(u)
          
          n=length(u);
          rplus=(u-circshift(u,1))./(circshift(u,-1)-u);
          rminus=1./rplus;
          
          minmodplus=max(0,min(1,rminus));
          minmodminus=max(0,min(1,rplus));
          
          ir=u+0.5*minmodplus.*(u-circshift(u,+1));   
          il=u-0.5*minmodminus.*(circshift(u,-1)-u);

 
          
      end
  
  function [ir, il]=MUSCLreconstructionvel(u,rhoir,rhoil,rho)
          
          n=length(u);
          rplus=(u-circshift(u,1))./(circshift(u,-1)-u);
          rminus=1./rplus;
          
          minmodplus=max(0,min(1,rminus));
          minmodminus=max(0,min(1,rplus));
          
          ir=u+0.5*(rhoil./rho).*minmodplus.*(u-circshift(u,+1));   
          il=u-0.5*(rhoir./rho).*minmodminus.*(circshift(u,-1)-u);
          
  end
  
  
 function [ur, ul]=MUSCLreconstruction2(u)
          
          n=length(u);
          
          
          du=minmod(u-circshift(u,1),circshift(u,-1)-u);
          
          ur=u+0.5*du;
          ul=u-0.5*du;
          
 end

 function [ur, ul]=MUSCLreconstruction2vel(u,rhoir,rhoil,rho)
          
          n=length(u);
         
          du=minmod(u-circshift(u,1),circshift(u,-1)-u);
          
          ur=u+0.5*du.*rhoil./rho;
          ul=u-0.5*du.*rhoir./rho;
          
 end
  


      function output=minmod(x,y) 
          
          output=zeros(length(x),1);
          for i=1:length(x)
              if x(i)>0 && y(i)>0
                  output(i)=min(x(i),y(i));
              elseif x(i)<0 && y(i)<0
                  output(i)=max(x(i),y(i));
              else
                  output(i)=0;
              end
          end
      end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % FUNCTIONS TO COMPUTE THE CONVOLUTIONS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
      function conv=kconvrhofunction(x,rho)
          deltax=x(2)-x(1);
          conv=deltax*sum((repmat(x',length(x),1)-x).^2/2.*rho);
          
%           conv=zeros(n,1);
%           for i=1:n/2
%               conv(i)=deltax*sum((x-x(i)).^2/2.*rho);
%           end
%           conv(n/2+1:end)=flipud(conv(1:n/2));
%           conv=conv';
      end
  
      function conv=kderconvrhofunction(x,rho)
          deltax=x(2)-x(1);
          conv=deltax*sum((repmat(x',length(x),1)-x).*rho);
      end
  
      function conv=CSconvfunction(x,rho,u)
          deltax=x(2)-x(1);
          conv=deltax*sum((1+abs(repmat(x',length(x),1)-x)).^(-0.25).*(repmat(u',length(x),1)-u(1:n)).*rho(1:n).*repmat(rho',length(x),1));
          
%           conv=zeros(n,1);
%           for i=1:n/2
%               conv(i)=-deltax*sum(rho(i)*(1+abs(x-x(i))).^(-0.25).*(u-u(i)).*rho);
%           end
%           conv(n/2+1:end)=-flipud(conv(1:n/2));
          %conv=conv';
      end
  end
  
  
