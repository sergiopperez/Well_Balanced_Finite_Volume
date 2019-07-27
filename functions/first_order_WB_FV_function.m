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
% FILE: first_order_WB_FV_function.m
%
% DESCRIPTION: compute the temporal derivative of the variable vector U
% with a first-oder well-balanced scheme, for a variety of free-energy
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


function dUdt = first_order_WB_FV_function(x,deltax,U,pd,nu,cik,alpha,cep,a,b,cefrHR,M,gamma,cCS,bc) 

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
% Construction of plus and minus values of the variable vector U at the 
% boundaries of the cells, by means of the well-balanced procedure 
%--------------------------------------------------------------------------
Uiplushalfplus=zeros(2*n,1);
Uiplushalfminus=zeros(2*n,1);

if pd>1    
    Uiplushalfplus(1:n)=((pd-1)/pd*max(((pd/(pd-1)*circshift(U(1:n),-1).^(pd-1)-max(H,circshift(H,-1))./nu+circshift(H,-1)./nu)),0)).^(1/(pd-1));
    Uiplushalfminus(1:n)=((pd-1)/pd*max(((pd/(pd-1)*U(1:n).^(pd-1)-max(H,circshift(H,-1))./nu+H./nu)),0)).^(1/(pd-1));
elseif pd==1
    Uiplushalfplus(1:n)=circshift(U(1:n),-1).*exp((-max(H,circshift(H,-1))+circshift(H,-1))./nu);
    Uiplushalfminus(1:n)=U(1:n).*exp((-max(H,circshift(H,-1))+H)./nu);
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
S=zeros(2*n,1);

% a) Part of the sources related to the external, interaction potentials
% and linear damping

if pd>1
    S(n+1:2*n)=1./deltax.*(nu*Uiplushalfminus(1:n).^pd-nu*circshift(Uiplushalfplus(1:n),+1).^pd)-gamma*U(1+n:2*n);
elseif pd==1
    S(n+1:2*n)=1./deltax.*(nu*Uiplushalfminus(1:n)-nu*circshift(Uiplushalfplus(1:n),+1))-gamma*U(1+n:2*n);    
end

% b) Part of the sources related to the Cucker-Smale damping term

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

    function conv=Wconvrhofunction(x,rho,alpha,deltax)   
        if alpha==0 % For convention W(x)=ln(x)
           conv=sum(log((repmat(x',length(x),1)-x)).*rho.*deltax); 
        elseif alpha>0
            conv=sum(abs((repmat(x',length(x),1)-x)).^alpha./alpha.*rho.*deltax)';
        elseif alpha<0 % Singularity treated by solving the integral analytically
            matrix1=abs(repmat(x',length(x),1)-x).^(alpha);
            matrix1(isinf(matrix1))=0;
            conv=sum(matrix1/(alpha).*rho.*deltax)'+(2/(alpha)/(alpha+1)*(deltax/2).^(alpha+1).*rho);
        end
    end

    function conv=Wgaussianconvrhofunction(x,rho,deltax)
        conv=sum(-exp(-abs((repmat(x',length(x),1)-x).^2./2)).*rho.*deltax./(2*pi))';
    end

    function conv=CSconvfunction(x,rho,rhou,deltax) 
        u=rhou./rho;
        conv=sum((1+abs(repmat(x',length(x),1)-x)).^(-0.25).*(repmat(u',length(x),1)-u(1:n)).*rho(1:n).*repmat(rho',length(x),1).*deltax)';
    end 

end