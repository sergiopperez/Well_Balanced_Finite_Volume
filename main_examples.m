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
% FILE: examples_main.m
%
% DESCRIPTION: compute the temporal derivative of the variable vector U
% with a first-oder well-balanced scheme, for a variety of free-energy
% funcitonal choices.


choiceorder=1; 
numberexample=3;

savefilm=0;
savescreenshots=0;
saveplots=1;

addpath('./functions/');


%--------------------------------------------------------------------------
% EXAMPLE 1: IDEAL-GAS PRESSURE AND ATTRACTIVE POTENTIAL
%--------------------------------------------------------------------------
if numberexample==1
n=200;
xboundary=linspace(-5,5,n+1)';
U0=zeros(2*n,1);
for j=1:n 
    U0(j)=integral(@(y)(0.2+cos(pi*(y)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)0.2+cos(pi*(y)/(xboundary(end)-xboundary(1))),xboundary(1),xboundary(end),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));  
    U0(j+n)=integral(@(y)(-0.05*sin(pi*(y)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
end
 pd=1;nu=1;cik=0;alpha=1;cep=1;a=0;b=0.5;gamma=1;cefrHR=0;cCS=0;timemax=15;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)
 
 
%--------------------------------------------------------------------------
% EXAMPLE 2: IDEAL-GAS PRESSURE, ATTRACTIVE POTENTIAL AND CUCKER_SMALE DAMPING TERM
%--------------------------------------------------------------------------
elseif numberexample==2
n=200;
xboundary=linspace(-5,5,n+1)';
U0=zeros(2*n,1);
for j=1:n 
    U0(j)=integral(@(y)(0.2+cos(pi*(y)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)0.2+cos(pi*(y)/(xboundary(end)-xboundary(1))),xboundary(1),xboundary(end),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));  
    U0(j+n)=integral(@(y)(-0.05*sin(pi*(y)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
end

 pd=1;nu=1;cik=0;alpha=1;cep=1;a=0;b=0.5;gamma=0;cefrHR=0;cCS=1;timemax=15;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)
 
%--------------------------------------------------------------------------
% EXAMPLE 3: IDEAL-GAS PRESSURE AND ATTRACTIVE KERNEL
%--------------------------------------------------------------------------
elseif numberexample==3
n=200;
xboundary=linspace(-5,5,n+1)';
U0=zeros(2*n,1);
for j=1:n 
    %U0(j)=integral(@(y)(0.2+cos(pi*(y)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)0.2+cos(pi*(y)/(xboundary(end)-xboundary(1))),xboundary(1),xboundary(end),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));  
    U0(j)=integral((@(y)exp(-y.^2/2)),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral((@(y)exp(-y.^2/2)),xboundary(1),xboundary(end),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));        
    %U0(j+n)=integral(@(y)(-0.05*sin(pi*(y)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
end
 pd=1;nu=1;cik=1;alpha=2;cep=0;a=0;b=0;gamma=0.01;cefrHR=0;cCS=0;timemax=75;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)

%--------------------------------------------------------------------------
% EXAMPLE 4: PRESSURE PROPORTIONAL TO THE SQUARE OF THE DENSITY AND
% ATTRACTIVE POTENTIAL
%--------------------------------------------------------------------------
elseif numberexample==4
n=200;
xboundary=linspace(-5,5,n+1)';
U0=zeros(2*n,1);
for j=1:n 
    U0(j)=integral(@(y)0.1+exp(-(y).^2),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)0.1+exp(-(y).^2),xboundary(1),xboundary(end),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
    U0(j+n)=integral(@(y)(-0.2*sin(pi*(y)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
end
 pd=2;nu=1;cik=0;alpha=2;cep=1;a=0;b=0.5;gamma=1;cefrHR=0;cCS=0;timemax=15;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)

%--------------------------------------------------------------------------
% EXAMPLE 5: MOVING WATER
%--------------------------------------------------------------------------
elseif numberexample==5
n=200;
xboundary=linspace(-8,9,n+1)';
x=(xboundary(2:end)+xboundary(1:end-1))/2;
U0=zeros(2*n,1);
for j=1:n 
    U0(j)=integral(@(y)(exp(-y.^2/2)),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)exp(-y.^2/2),xboundary(1),xboundary(end))/(xboundary(j+1)-xboundary(j));
    U0(j+n)=U0(j)*0.2;
end
 pd=1;nu=1;cik=1;alpha=2;cep=0;a=0;b=0;gamma=0;cefrHR=0;cCS=1;timemax=3;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)
 

%--------------------------------------------------------------------------
% EXAMPLE 6.1: PRESSURE PROPORTIONAL TO THE SQUARE OF THE DENSITY AND
% DOUBLE-WELL POTENTIAL. SYMMETRIC STEADY DENSITY
%--------------------------------------------------------------------------
elseif numberexample==61
n=200;
xboundary=linspace(-5,5,n+1)';
U0=zeros(2*n,1);
xmean=0;
for j=1:n 
    U0(j)=integral(@(y)0.1+exp(-(y-xmean).^2),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)0.1+exp(-(y-xmean).^2),xboundary(1),xboundary(end),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
    U0(j+n)=integral(@(y)(-0.2*sin(pi*(y)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
end
 pd=2;nu=1;cik=0;alpha=2;cep=1;a=0.25;b=-1.5;gamma=1;cefrHR=0;cCS=0;timemax=20;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)

%--------------------------------------------------------------------------
% EXAMPLE 6.2: PRESSURE PROPORTIONAL TO THE SQUARE OF THE DENSITY AND
% DOUBLE-WELL POTENTIAL. ASYMMETRIC STEADY DENSITY
%--------------------------------------------------------------------------
elseif numberexample==62
n=200;
xboundary=linspace(-5,5,n+1)';
U0=zeros(2*n,1);
xmean=1;
for j=1:n 
    U0(j)=integral(@(y)0.1+exp(-(y-xmean).^2),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)0.1+exp(-(y-xmean).^2),xboundary(1),xboundary(end),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
    U0(j+n)=integral(@(y)(-0.2*sin(pi*(y)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
end
 pd=2;nu=1;cik=0;alpha=2;cep=1;a=0.25;b=-1.5;gamma=1;cefrHR=0;cCS=0;timemax=20;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)
 
%--------------------------------------------------------------------------
% EXAMPLE 6.3: PRESSURE PROPORTIONAL TO THE SQUARE OF THE DENSITY AND
% DOUBLE-WELL POTENTIAL. SYMMETRIC STEADY DENSITY (2)
%--------------------------------------------------------------------------
elseif numberexample==63
n=200;
xboundary=linspace(-5,5,n+1)';
U0=zeros(2*n,1);
xmean=1;
for j=1:n 
    U0(j)=integral(@(y)0.1+exp(-(y-xmean).^2),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)0.1+exp(-(y-xmean).^2),xboundary(1),xboundary(end),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
    U0(j+n)=integral(@(y)(-0.2*sin(pi*(y)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));
end
 pd=2;nu=1;cik=0;alpha=2;cep=1;a=0.25;b=-0.5;gamma=1;cefrHR=0;cCS=0;timemax=20;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)

Un nombre très limité de verbes au futur simple et au futur proche. (Est-ce parce que la tentative d’accéder à une vie différente, voire meilleure, des clandestins est vouée à l’échec ?)

%--------------------------------------------------------------------------
% EXAMPLE 7: IDEAL PRESSURE WITH NOISE PARAMETER AND ITS PHASE TRANSITION
%--------------------------------------------------------------------------
elseif numberexample==7
n=200;
xboundary=linspace(-5,5,n+1)';
U0=zeros(2*n,1);
xmean=0.1;
for j=1:n 
    U0(j)=integral(@(y)(0.2+cos(pi*(y-xmean)/(xboundary(end)-xboundary(1)))),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)0.2+cos(pi*(y-xmean)/(xboundary(end)-xboundary(1))),xboundary(1),xboundary(end),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));  
    U0(j+n)=0;
end
 pd=1;nu=0.4;cik=1;alpha=2;cep=1;a=0.25;b=-0.5;gamma=3;cefrHR=0;cCS=0;timemax=13;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)

%--------------------------------------------------------------------------
% EXAMPLE 8.1: KELLER-SEGEL SYSTEM: COMPACTLY-SUPPORTED STEADY STATE
%--------------------------------------------------------------------------
elseif numberexample==81
n=200;
xboundary=linspace(-8,8,n+1)';
U0=zeros(2*n,1);
for j=1:n 
    U0(j)=integral(@(y)(exp(-4*(y+2).^2/10)+exp(-4*(y-2).^2/10)),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)(exp(-4*(y+2).^2/10)+exp(-4*(y-2).^2/10)),-1000,1000,'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));   
    U0(j+n)=0;
end
 pd=1.5;nu=1;cik=1;alpha=0.5;cep=0;a=0;b=0;gamma=1;cefrHR=0;cCS=0;timemax=70;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)


 
%--------------------------------------------------------------------------
% EXAMPLE 8.2: KELLER-SEGEL SYSTEM: FINITE-TIME BLOW UP
%--------------------------------------------------------------------------
elseif numberexample==82
n=201;
xboundary=linspace(-8,8,n+1)';
U0=zeros(2*n,1);
for j=1:n 
    U0(j)=integral(@(y)(exp(-4*(y+2).^2/10)+exp(-4*(y-2).^2/10)),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)(exp(-4*(y+2).^2/10)+exp(-4*(y-2).^2/10)),-1000,1000,'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));   
    U0(j+n)=0;
end
 pd=1.3;nu=1;cik=1;alpha=-0.5;cep=0;a=0;b=0;gamma=1;cefrHR=0;cCS=0;timemax=250;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)

%--------------------------------------------------------------------------
% EXAMPLE 8.3: KELLER-SEGEL SYSTEM: MORSE-TYPE POTENTIAL
%--------------------------------------------------------------------------
elseif numberexample==83
n=200;
xboundary=linspace(-8,12,n+1)';
U0=zeros(2*n,1);
for j=1:n 
    U0(j)=1.2*integral(@(y)(exp(-0.5*(y+3).^2)+exp(-0.5*(y-3).^2)+0.55*exp(-0.5*(y-8.5).^2)),xboundary(j),xboundary(j+1),'RelTol',1e-30)/integral(@(y)(exp(-0.5*(y+3).^2)+exp(-0.5*(y-3).^2)+0.55*exp(-0.5*(y-8.5).^2)),-1000,1000,'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));   
    U0(j+n)=0;
end
 pd=3;nu=1;cik=2;alpha=0;cep=0;a=0;b=0;gamma=0.05;cefrHR=0;cCS=0;timemax=3000;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)

%--------------------------------------------------------------------------
% EXAMPLE 9.1: DDFT HARD RODS WITH CONFINING POTENTIAL
%--------------------------------------------------------------------------
elseif numberexample==91
n=200;
xboundary=[linspace(-13,-6.01,30) linspace(-6,6,141) linspace(6.01,13,30)]';     
U0=zeros(2*n,1);
for j=1:n 
    U0(j)=integral(@(y)(exp(-y.^2/20.372)),xboundary(j),xboundary(j+1),'RelTol',1e-30)/(xboundary(j+1)-xboundary(j));   
    U0(j+n)=0;bbjhcjdjbjvjvnm    nvknkvbnvknvknh####5363
end
 pd=1;nu=1;cik=0;alpha=0;cep=1;a=0;b=1;gamma=1;cefrHR=1;cCS=0;timemax=20;bc=1;
 RK3TVD_timeinteg_function(n,xboundary,U0,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)

 %--------------------------------------------------------------------------
% EXAMPLE 9.2: DDFT HARD RODS AND DIFFUSION
%--------------------------------------------------------------------------
elseif numberexample==92
n=200;
xboundary=[linspace(-13,-6.01,30) linspace(-6,6,141) linspace(6.01,13,30)]';     
load([pwd '/data/IC-92.mat'])
pd=1;nu=1;cik=0;alpha=0;cep=1;a=0;b=0;gamma=1;cefrHR=1;cCS=0;timemax=30;bc=1;
RK3TVD_timeinteg_function(n,xboundary,Usteady81,pd,nu,cik,alpha,cep,a,b,gamma,cefrHR,cCS,timemax,choiceorder,bc,numberexample,savefilm,savescreenshots,saveplots)


end 