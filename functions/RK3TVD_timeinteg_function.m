%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CODE FOR THE PAPER "WELL-BALANCED FINITE VOLUME SCHEMES FOR HYDRODYNAMIC
% EQUATIONS WITH GENERAL FREE ENERGY"
%
% 1st AND 2nd ORDER 1D WELL-BALANCED FINITE VOLUMES
% 
% AUTHOR OF THE CODE: SERGIO P. PEREZ
%
% COAUTHORS: JOSÉ A. CARRILLO, SERAFIM KALLIADASIS, CHI-WANG SHU
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FILE: RK3TVD_timeinteg_function.m
%
% DESCRIPTION: compute the temporal evolution by means of a total variation
% diminishing Runge-Kutta of third order, inspired in the work of 
% S. Gottlieb and C.-W. Shu, Total variation diminishing Runge-Kutta schemes, 
% Math. Comput. Am. Math. Soc., 67 (1998), pp. 73--85.
%
% INPUTS:
%     n-> number of cells in the domain
%     xboundary-> vector of length n+1 containing the positions of the cell
%         interfaces
%     U0-> vector of length 2*n with the initial conditions for density and
%         momentum
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
%     gamma-> coefficient of the linear damping
%     cefrHR-> logical value representing the choice of including the
%              excessive free energy of hard rods
%     cCS-> logical value representing the choice of including the
%         Cucker-Smale damping terms from collective behaviour
%     timemax-> time at which the simulation finishes
%     choiceorder-> select between first order (choiceorder==1) or second
%         order (choiceorder==2)
%     bc-> value indicating the chosen boundary conditions. If bc==1, the
%         boundary conditions are no flux, and if bc==2 they are periodic
%     numberexample-> number of the example from the paper
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)

% Create vector of positions of the central node at each cell
x=(xboundary(2:end)+xboundary(1:end-1))/2;

% Create vector with the length of each cell
deltax=xboundary(2:end)-xboundary(1:end-1);

% Select CFL condition
cfl=0.2;

% Create matrix employed for the rods
if cefrHR==1
    sigma = 1;
    M = Fplus(x,sigma);
else
    M=1;
end

% Start simulation

ntimes=100000;
U=zeros(2*n,ntimes);
U(:,1)=U0;
t=zeros(ntimes,1);
       
i=1;
while t(i)<timemax

    % Select time step 
    if pd>1
         deltat=cfl*min(deltax)/max(abs(U(n+1:2*n,i)./U(1:n,i))+3^((pd-1)/4));      
    elseif pd==1      
        lambda=max(abs(U(n+1:2*n,i)./U(1:n,i)+sqrt(nu)),abs(U(n+1:2*n,i)./U(1:n,i)-sqrt(nu)));       
        maxlambda=max(lambda);  
        deltat=cfl*min(deltax)/maxlambda;       
    end
    
    % Advance in time
    t(i+1)=t(i)+deltat;
    
    % Apply the RK3 TVD
    
    if choiceorder==1      
        psi_1=U(:,i)+deltat*first_order_WB_FV_function(x,deltax,U(:,i),pd,nu,cik,alpha,cep,a,b,cefrHR,M,gamma,cCS,bc);       
        psi_2=3/4*U(:,i)+1/4*psi_1+1/4*deltat*first_order_WB_FV_function(x,deltax,psi_1,pd,nu,cik,alpha,cep,a,b,cefrHR,M,gamma,cCS,bc);       
        U(:,i+1)=1/3*U(:,i)+2/3*psi_2+2/3*deltat*first_order_WB_FV_function(x,deltax,psi_2,pd,nu,cik,alpha,cep,a,b,cefrHR,M,gamma,cCS,bc);
        
    elseif choiceorder==2      
        psi_1=U(:,i)+deltat*second_order_FV_LLF_FUNCTION(U(:,i),deltax,x',choicewellbalance,choicepressure,choicekernel,choicepotential,choiceCS,V);     
        psi_2=3/4*U(:,i)+1/4*psi_1+1/4*deltat*second_order_FV_LLF_FUNCTION(psi_1,deltax,x',choicewellbalance,choicepressure,choicekernel,choicepotential,choiceCS,V);      
        U(:,i+1)=1/3*U(:,i)+2/3*psi_2+2/3*deltat*second_order_FV_LLF_FUNCTION(psi_2,deltax,x',choicewellbalance,choicepressure,choicekernel,choicepotential,choiceCS,V);       
    end
    
    % Message in command window
    display(['--------------------'])
    display(['Time: ',num2str(t(i))])
    display(['L1 norm of the difference between the new and old state: ',num2str(norm(U(:,i+1)-U(:,i),1))])
    
    % Advance to next time index
    i=i+1;
end

% Remove the columns that have not been filled
ntimes=i;
U=U(:,1:ntimes);
t=t(1:ntimes);

% Call the function producefigures
if saveplots==1
    producefigures(n,x,deltax,U,t,pd,nu,alpha,a,b,M,numberexample)
end

% Call the function producefilms
if savefilm==1 || savescreenshots==1
    producefilms(n,x,U,t,numberexample)
end

end