function producefigures(n,x,deltax,U,t,pd,nu,alpha,a,b,M,numberexample)

close all

%--------------------------------------------------------------------------
% EXAMPLE 1: ideal-gas pressure and ideal potential
%--------------------------------------------------------------------------
if numberexample==1
    times=['$$t=0$$  '; '$$t=0.7$$' ;'$$t=2$$  '; '$$t=15$$ '];
    [fff,index1]=min(abs(t(:)-0.7));[fff,index2]=min(abs(t(:)-2));[fff,index3]=min(abs(t(:)-15));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=nu*log(U(1:n,timespos(cases)))+b*x.^2;
    end
    freeenergy=sum(deltax.*(nu*U(1:n,:).*(log(U(1:n,:))-1)+b*U(1:n,:).*x.^2));
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
    
%--------------------------------------------------------------------------
% EXAMPLE 2: ideal-gas pressure, attractive potential and Cucker-Smale
% damping term
%--------------------------------------------------------------------------
elseif numberexample==2;
    times=['$$t=0$$  '; '$$t=0.7$$' ;'$$t=2$$  '; '$$t=15$$ '];
    [fff,index1]=min(abs(t(:)-0.7));[fff,index2]=min(abs(t(:)-2));[fff,index3]=min(abs(t(:)-15));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=nu*log(U(1:n,timespos(cases)))+b*x.^2;
    end
    freeenergy=sum(deltax.*(nu*U(1:n,:).*(log(U(1:n,:))-1)+b*U(1:n,:).*x.^2));
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
    
%--------------------------------------------------------------------------
% EXAMPLE 3: ideal-gas pressure and attractive potential
%--------------------------------------------------------------------------
elseif numberexample==3;
    times=['$$t=0$$  '; '$$t=2.5$$' ;'$$t=5$$  '; '$$t=75$$ '];
    [fff,index1]=min(abs(t(:)-2.5));[fff,index2]=min(abs(t(:)-5));[fff,index3]=min(abs(t(:)-75));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=nu*log(U(1:n,timespos(cases)))+sum(abs(repmat(x',length(x),1)-x).^alpha./alpha.*U(1:n,timespos(cases)).*deltax)';
    end
    conv=zeros(length(t),1);
    for nt=1:length(t)
        conv(nt)=0.5*sum(deltax.*U(1:n,nt).*sum(abs(repmat(x',length(x),1)-x).^alpha./alpha.*U(1:n,nt).*deltax)');
    end
    freeenergy=sum(deltax.*(nu*U(1:n,:).*(log(U(1:n,:))-1)))+conv';
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
    
%--------------------------------------------------------------------------
% EXAMPLE 4: pressure porportional to the square of the density and 
% attractive potential
%--------------------------------------------------------------------------
elseif numberexample==4;
    times=['$$t=0$$  '; '$$t=0.7$$' ;'$$t=3$$  '; '$$t=15$$ '];
    [fff,index1]=min(abs(t(:)-0.7));[fff,index2]=min(abs(t(:)-3));[fff,index3]=min(abs(t(:)-15));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=pd/(pd-1)*U(1:n,timespos(cases))+b*x.^2;
    end
    freeenergy=sum(deltax.*(nu/(pd-1)*U(1:n,:).^(pd)+b*U(1:n,:).*x.^2));
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)

%--------------------------------------------------------------------------
% EXAMPLE 5: MOVING WATER
%--------------------------------------------------------------------------
elseif numberexample==5;
    times=['$$t=0$$  '; '$$t=1$$  ' ;'$$t=2$$  '; '$$t=3$$  '];
    [fff,index1]=min(abs(t(:)-1));[fff,index2]=min(abs(t(:)-2));[fff,index3]=min(abs(t(:)-3));
    timespos=[1;index1;index2;index3];    

        variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=nu*log(U(1:n,timespos(cases)))+sum(abs(repmat(x',length(x),1)-x).^alpha./alpha.*U(1:n,timespos(cases)).*deltax)';
    end
    conv=zeros(length(t),1);
    for nt=1:length(t)
        conv(nt)=0.5*sum(deltax.*U(1:n,nt).*sum(abs(repmat(x',length(x),1)-x).^alpha./alpha.*U(1:n,nt).*deltax)');
    end
    freeenergy=sum(deltax.*(nu*U(1:n,:).*(log(U(1:n,:))-1)))+conv';
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
    

%--------------------------------------------------------------------------
% EXAMPLE 6.1: pressure porportional to the square of the density and 
% double-well potential, symmetric steady density
%--------------------------------------------------------------------------
elseif numberexample==61;
    times=['$$t=0$$  '; '$$t=0.3$$' ;'$$t=2.5$$'; '$$t=20$$ '];
    [fff,index1]=min(abs(t(:)-0.3));[fff,index2]=min(abs(t(:)-2.5));[fff,index3]=min(abs(t(:)-20));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=pd/(pd-1)*U(1:n,timespos(cases))+a*x.^4+b*x.^2;
    end
    freeenergy=sum(deltax.*(nu/(pd-1)*U(1:n,:).^(pd)+a*U(1:n,:).*x.^4+b*U(1:n,:).*x.^2));
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
    
%--------------------------------------------------------------------------
% EXAMPLE 6.2: pressure porportional to the square of the density and 
% double-well potential, asymmetric steady state
%--------------------------------------------------------------------------
elseif numberexample==62;
    times=['$$t=0$$  '; '$$t=0.3$$' ;'$$t=2.5$$'; '$$t=20$$ '];
    [fff,index1]=min(abs(t(:)-0.3));[fff,index2]=min(abs(t(:)-2.5));[fff,index3]=min(abs(t(:)-20));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=pd/(pd-1)*U(1:n,timespos(cases))+a*x.^4+b*x.^2;
    end
    freeenergy=sum(deltax.*(nu/(pd-1)*U(1:n,:).^(pd)+a*U(1:n,:).*x.^4+b*U(1:n,:).*x.^2));
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
    
%--------------------------------------------------------------------------
% EXAMPLE 6.3: pressure porportional to the square of the density and 
% double-well potential, symmetric steady state(2)
%--------------------------------------------------------------------------
elseif numberexample==63;
    times=['$$t=0$$  '; '$$t=0.3$$' ;'$$t=2.5$$'; '$$t=20$$ '];
    [fff,index1]=min(abs(t(:)-0.3));[fff,index2]=min(abs(t(:)-2.5));[fff,index3]=min(abs(t(:)-20));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=pd/(pd-1)*U(1:n,timespos(cases))+a*x.^4+b*x.^2;
    end
    freeenergy=sum(deltax.*(nu/(pd-1)*U(1:n,:).^(pd)+a*U(1:n,:).*x.^4+b*U(1:n,:).*x.^2));
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)

%--------------------------------------------------------------------------
% EXAMPLE 7: pressure porportional to the square of the density and 
% double-well potential, symmetric steady density
%--------------------------------------------------------------------------
elseif numberexample==7;
    times=['$$t=0$$  '; '$$t=0.3$$' ;'$$t=2.5$$'; '$$t=20$$ '];
    [fff,index1]=min(abs(t(:)-0.3));[fff,index2]=min(abs(t(:)-2.5));[fff,index3]=min(abs(t(:)-50));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=nu*log(U(1:n,timespos(cases)))+sum(abs(repmat(x',length(x),1)-x).^alpha./alpha.*U(1:n,timespos(cases)).*deltax)'+a*x.^4+b*x.^2;
    end
    conv=zeros(length(t),1);
    for nt=1:length(t)
        conv(nt)=0.5*sum(deltax.*U(1:n,nt).*sum(abs(repmat(x',length(x),1)-x).^alpha./alpha.*U(1:n,nt).*deltax)');
    end
    freeenergy=sum(deltax.*(nu*U(1:n,:).*(log(U(1:n,:))-1)+a*U(1:n,:).*x.^4+b*U(1:n,:).*x.^2))+conv';
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
%--------------------------------------------------------------------------
% EXAMPLE 8.1: Keller-Segel system: compactly-supported steady state
%--------------------------------------------------------------------------
elseif numberexample==81;
    
    times=['$$t=0$$  '; '$$t=3.7$$' ;'$$t=11$$ '; '$$t=70$$ '];
    [fff,index1]=min(abs(t(:)-3.7));[fff,index2]=min(abs(t(:)-11));[fff,index3]=min(abs(t(:)-70));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=pd/(pd-1)*U(1:n,timespos(cases)).^(pd-1)+sum(abs(repmat(x',length(x),1)-x).^alpha./alpha.*U(1:n,timespos(cases)).*deltax)';
    end
    conv=zeros(length(t),1);
    for nt=1:length(t)
        conv(nt)=0.5*sum(deltax.*U(1:n,nt).*sum(abs(repmat(x',length(x),1)-x).^alpha./alpha.*U(1:n,nt).*deltax)');
    end
    
    freeenergy=sum(deltax.*(nu/(pd-1)*U(1:n,:).^(pd)))+conv';
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
    
%--------------------------------------------------------------------------
% EXAMPLE 8.2: Keller-Segel system: finite-time blow up
%--------------------------------------------------------------------------
elseif numberexample==82;
    
    times=['$$t=0$$  '; '$$t=9$$  ' ;'$$t=65$$ '; '$$t=250$$'];
    [fff,index1]=min(abs(t(:)-9));[fff,index2]=min(abs(t(:)-65));[fff,index3]=min(abs(t(:)-900));
    timespos=[1;index1;index2;index3];  
    
    variationfreeenergy=zeros(n,4);
    matrix1=abs(repmat(x',length(x),1)-x).^(alpha)./alpha; 
    matrix1(isinf(matrix1))=0;
    for cases=1:4
        conv=sum(matrix1.*U(1:n,timespos(cases)).*deltax)'+(2/(alpha)/(alpha+1)*(deltax/2).^(alpha+1).*U(1:n,timespos(cases)));
        variationfreeenergy(1:n,cases)=pd/(pd-1)*U(1:n,timespos(cases)).^(pd-1)+conv;
    end
    conv=zeros(length(t),1);
    for nt=1:length(t)
        conv(nt)=0.5*sum(deltax.*U(1:n,nt).*(sum(matrix1.*U(1:n,timespos(cases)).*deltax)'+(2/(alpha)/(alpha+1)*(deltax/2).^(alpha+1).*U(1:n,timespos(cases)))));
    end 
    freeenergy=sum(deltax.*(nu/(pd-1)*U(1:n,:).^(pd)))+conv';
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
    createfigurealone(x,U(1:n,timespos(4)),times,numberexample)
    
%--------------------------------------------------------------------------
% EXAMPLE 8.3: pressure proportional to the cube of the density and Morse-type potential 
%--------------------------------------------------------------------------
elseif numberexample==83;
    times=['$$t=0$$   '; '$$t=100$$ ' ;'$$t=270$$ '; '$$t=3000$$'];
    [fff,index1]=min(abs(t(:)-100));[fff,index2]=min(abs(t(:)-270));[fff,index3]=min(abs(t(:)-3000));
    timespos=[1;index1;index2;index3];
     
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=nu*pd/(pd-1)*U(1:n,timespos(cases)).^(pd-1)+sum(-exp(-abs((repmat(x',length(x),1)-x).^2./2)).*U(1:n,timespos(cases)).*deltax./(2*pi))';
    end
    conv=zeros(length(t),1);
    for nt=1:length(t)
        conv(nt)=0.5*sum(deltax.*U(1:n,nt).*sum(-exp(-abs((repmat(x',length(x),1)-x).^2./2)).*U(1:n,nt).*deltax./(2*pi))');
    end
    freeenergy=sum(deltax.*(nu/(pd-1)*U(1:n,:).^(pd)))+conv';
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
    [fff,index4]=min(abs(t(:)-400))
    createfigure3D(x,U,t,numberexample,index4)
    
%--------------------------------------------------------------------------
% EXAMPLE 9.1: DDFT hard rods with confining potential
%--------------------------------------------------------------------------
    
elseif numberexample==91
    times=['$$t=0$$  '; '$$t=0.6$$' ;'$$t=2.3$$'; '$$t=20$$ '];
    [fff,index1]=min(abs(t(:)-0.5));[fff,index2]=min(abs(t(:)-2.3));[fff,index3]=min(abs(t(:)-20));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=nu*log(U(1:n,timespos(cases)))-log(1-M*U(1:n,timespos(cases)))+rot90(M,2)*(U(1:n,timespos(cases))./(1-M*U(1:n,timespos(cases))))+b*(x.^2);
        
    end
    freeenergy=sum(deltax.*(nu*U(1:n,:).*(log(U(1:n,:))-1)+b*U(1:n,:).*x.^2-0.5*U(1:n,:).*log(1-M*U(1:n,:))-0.5*U(1:n,:).*log(1-rot90(M,2)*U(1:n,:))));
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)
    
    % Save last state for example 8.2
    Usteady91=U(:,end);
    save([pwd '/data/IC-92.mat'],'Usteady91');

%--------------------------------------------------------------------------
% EXAMPLE 9.2: DDFT hard rods diffusing without confining potential
%--------------------------------------------------------------------------
    
elseif numberexample==92
    times=['$$t=0$$  '; '$$t=0.6$$' ;'$$t=2.3$$'; '$$t=30$$ '];
    [fff,index1]=min(abs(t(:)-0.5));[fff,index2]=min(abs(t(:)-2.3));[fff,index3]=min(abs(t(:)-30));
    timespos=[1;index1;index2;index3];
    
    variationfreeenergy=zeros(n,4);
    for cases=1:4
        variationfreeenergy(:,cases)=nu*log(U(1:n,timespos(cases)))-log(1-M*U(1:n,timespos(cases)))+rot90(M,2)*(U(1:n,timespos(cases))./(1-M*U(1:n,timespos(cases))))+b*(x.^2);
        
    end
    freeenergy=sum(deltax.*(nu*U(1:n,:).*(log(U(1:n,:))-1)+b*U(1:n,:).*x.^2-0.5*U(1:n,:).*log(1-M*U(1:n,:))-0.5*U(1:n,:).*log(1-rot90(M,2)*U(1:n,:))));
    totalenergy=sum(0.5*deltax.*U(n+1:2*n,:).^2./U(1:n,:))+freeenergy;
    
    createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,timespos(1),timespos(2),timespos(3),timespos(4),numberexample)   
    
end

