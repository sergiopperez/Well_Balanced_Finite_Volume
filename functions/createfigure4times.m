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
% FILE: createfigure4times.m
%
% DESCRIPTION: create four figures: density profile at four different
% times, momentum profile at four different times, variation of the free
% energy profile at four different times, and temporal evolution of the
% total and free energy
%
% INPUTS:
%     x-> 1D mesh with n nodes located at the centre of the finite volume 
%         cells
%     U-> matrix whose columns are a vector of variables (density and 
%         momentum) with length 2*n. Each column corresponds to a certain
%         time of the simulation
%     variationfreeenergy-> matrix with four columns of length n, with each
%         column containning the values of the variation of the free energy
%         with respect to the density across the domain. Each column refers
%         to a particular time in the simulation (t1, t2, t3 and t4)
%     freeenergy-> vector containing the values of the free energy at each
%         time of the simulation
%     totalenergy-> vector containing the values of the total energy at each
%         time of the simulation
%     t-> vector of times
%     times-> array of characters with four elements, each of them contains
%         the characters to indicate the time of each plot in the legend
%     t1-> index of time of the first plot
%     t2-> index of time of the second plot
%     t3-> index of time of the third plot
%     t4-> index of time of the forth plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

function createfigure4times(x,U,variationfreeenergy,freeenergy,totalenergy,t,times,t1,t2,t3,t4,number)

%--------------------------------------------------------------------------
% First plot: density profile at four different times
%-------------------------------------------------------------------------- 

% Extract the density at four different times from the matrix U
Y1=U(1:length(x),t1);
Y2=U(1:length(x),t2);
Y3=U(1:length(x),t3);
Y4=U(1:length(x),t4);

% Create figure
figure1 = figure('PaperOrientation','landscape','pos',[10 10 560 450],'PaperSize',[20 13]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Customize axes for each example
if number==1
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.6]);
    yticks([0 0.3 0.6]);
elseif number==2
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.9]);
    yticks([0 0.3 0.6 0.9]);
elseif number==3
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.6]);
    yticks([0 0.3 0.6]);
elseif number==4
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.6]);
    yticks([0 0.3 0.6]);
elseif number==5
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.5]);
    yticks([0 0.2 0.4]);  
elseif number==61 || number==62 || number==63
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.8]);
    yticks([0 0.4 0.8]);
elseif number==7
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.9]);
    yticks([0 0.3 0.6 0.9]);
elseif number==81
    xlim(axes1,[-8 8]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.4]);
    yticks([0 0.2 0.4]);
elseif number==82
    xlim(axes1,[-8 8]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.4]);
    yticks([0 0.2 0.4]);
elseif number==83
    xlim(axes1,[-8 12]);
    xticks([-6 0 6 12]);
    ylim(axes1,[0 0.24]);
    yticks([0 0.1 0.2]);
elseif number==91
    xlim(axes1,[-10 10]);
    xticks([-10 -5 0 5 10]);
    ylim(axes1,[0 2]);
    yticks([0 1 2]);
elseif number==92
    xlim(axes1,[-10 10]);
    xticks([-10 -5 0 5 10]);
    ylim(axes1,[0 2]);
    yticks([0 1 2]);  
end

% Create 4 plots
plot(x,Y1,'DisplayName',times(1,:),'LineWidth',3,'LineStyle',':',...
    'Color',[0 0.5647 0.6196]);

hold all
plot(x,Y2,'DisplayName',times(2,:),'LineWidth',3,'LineStyle','--',...
    'Color',[0.4980 0.6902 0.0196]);

hold all
plot(x,Y3,'DisplayName',times(3,:),'LineWidth',3,'LineStyle','-.',...
    'Color',[0.87058824300766 0.490196079015732 0]);

hold all
plot(x,Y4,'DisplayName',times(4,:),'LineWidth',3,'LineStyle','-',...
    'Color',[0.4863 0.0784 0.3020]);

% Set the remaining axes properties
box(axes1,'on');
set(axes1,'FontSize',16,'TickLabelInterpreter','latex','XGrid','on','YGrid',...
    'on');

% Create xlabel
xlabel({'$$x$$'},'FontSize',23,'Interpreter','latex');

% Create legend
if number==62
    legend1 = legend(axes1,'show','Location','northwest');
else
    legend1 = legend(axes1,'show','Location','northeast');
end
set(legend1,'Interpreter','latex','FontSize',21);

% Create textarrow for the ylabel
annotation(figure1,'textarrow',[0.033 0.1],...
    [0.56 0.865459627329192],'TextEdgeColor','none',...
    'FontSize',23,...
    'Interpreter','latex',...
    'String',{'$$\rho$$'},...
    'HeadStyle','none',...
    'LineStyle','none');

% Save the figure
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
saveas(gcf, [pwd sprintf('/figures/density-%d',number)], 'pdf')%Save figure
%--------------------------------------------------------------------------
% Second plot: momentum profile at four different times
%-------------------------------------------------------------------------- 

% Extract the momentum at four different times from the matrix U
Y5=U(length(x)+1:2*length(x),t1);
Y6=U(length(x)+1:2*length(x),t2);
Y7=U(length(x)+1:2*length(x),t3);
Y8=U(length(x)+1:2*length(x),t4);

% Uncomment the following lines if you prefer to plot the velocity
% Y5=U(length(x)+1:2*length(x),t1)./U(1:length(x),t1);
% Y6=U(length(x)+1:2*length(x),t2)./U(1:length(x),t2);
% Y7=U(length(x)+1:2*length(x),t3)./U(1:length(x),t3);
% Y8=U(length(x)+1:2*length(x),t4)./U(1:length(x),t4);

% Create figure
figure2 = figure('PaperOrientation','landscape','pos',[10 10 560 450],'PaperSize',[20 13]);

% Create axes
axes2 = axes('Parent',figure2);
hold(axes2,'on');

% Customize axes for each example
if number==1
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.2 0.2]);
    yticks([-0.2 0 0.2]);
elseif number==2
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.2 0.2]);
    yticks([-0.2 0 0.2]);
elseif number==3
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.35 0.35]);
    yticks([-0.3 0 0.3]);
elseif number==4
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.42 0.42]);
    yticks([-0.4 0 0.4]);
elseif number==5
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[0 0.09]);
    yticks([0 0.08]);
%     ylim(axes2,[0 0.4]);
%     yticks([0 0.2 0.4]);
elseif number==61 || number==62 || number==63
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.4 0.4]);
    yticks([-0.4 0 0.4]);
elseif number==7
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.2 0.2]);
    yticks([-0.2 0 0.2]);
elseif number==81
    xlim(axes2,[-8 8]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.04 0.04]);
    yticks([-0.04 0 0.04]);
elseif number==82
    xlim(axes2,[-8 8]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.08 0.08]);
    yticks([-0.08 0 0.08]);
elseif number==83
    xlim(axes2,[-8 12]);
    xticks([-6 0 6 12]);
    ylim(axes2,[-0.005 0.005]);
    yticks([-0.005 0 0.005])
elseif number==91
    xlim(axes2,[-10 10]);
    xticks([-10 -5 0 5 10]);
    ylim(axes2,[-1 1]);
    yticks([-1 0 1]);
elseif number==92
    xlim(axes2,[-10 10]);
    xticks([-10 -5 0 5 10]);
    ylim(axes2,[-1.3 1.3]);
    yticks([-1 0 1]);
end

% Create 4 plots
plot(x,Y5,'DisplayName',times(1,:),'LineWidth',3,'LineStyle',':',...
    'Color',[0 0.5647 0.6196]);

hold all
plot(x,Y6,'DisplayName',times(2,:),'LineWidth',3,'LineStyle','--',...
    'Color',[0.4980 0.6902 0.0196]);

hold all
plot(x,Y7,'DisplayName',times(3,:),'LineWidth',3,'LineStyle','-.',...
    'Color',[0.87058824300766 0.490196079015732 0]);

hold all
plot(x,Y8,'DisplayName',times(4,:),'LineWidth',3,'LineStyle','-',...
    'Color',[0.4863 0.0784 0.3020]);

% Set the remaining axes properties
box(axes2,'on');
set(axes2,'FontSize',16,'TickLabelInterpreter','latex','XGrid','on','YGrid',...
    'on');

% Create xlabel
xlabel({'$$x$$'},'FontSize',23,'Interpreter','latex');

% Create legend
if number==92
    legend2 = legend(axes2,'show','Location','northwest');
else
    legend2 = legend(axes2,'show','Location','northeast');
end
set(legend2,'Interpreter','latex','FontSize',21);

% Create textarrow for the ylabel
annotation(figure2,'textarrow',[0.033 0.1],...
    [0.56 0.865459627329192],'TextEdgeColor','none',...
    'FontSize',23,...
    'String',{'$$\rho u$$'},...
    'HeadStyle','none',...
    'LineStyle','none','Interpreter','latex');

% Save the figure
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
saveas(gcf, [pwd sprintf('/figures/momentum-%d',number)], 'pdf')%Save figure
%saveas(gcf, [pwd sprintf('/figures/velocity-%d',number)], 'pdf')%Save figure

%--------------------------------------------------------------------------
% Third plot: variation of the free energy profile at four different times
%-------------------------------------------------------------------------- 

% Create figure
figure3 = figure('PaperOrientation','landscape','pos',[10 10 560 450],'PaperSize',[20 13]);

% Create axes
axes3 = axes('Parent',figure3);
hold(axes3,'on');

% Customize axes for each example
if number==1
    xlim(axes3,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes3,[-4 10]);
    yticks([-4 0 4 8]);
elseif number==2
    xlim(axes3,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes3,[-4 10]);
    yticks([-4 0 4 8]);
elseif number==3
    xlim(axes3,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes3,[-10 13]);
    yticks([-10 0 10]);
elseif number==4
    xlim(axes3,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes3,[0 13]);
    yticks([0 6 12]);
elseif number==61 || number==62 || number==63
    xlim(axes3,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes3,[-3 8]);
    yticks([-3 0 3 6]);
elseif number==7
    xlim(axes3,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes3,[-4 10]);
    yticks([-4 0 4 8]);
elseif number==81
    xlim(axes3,[-8 8]);
    xticks([-5 0 5]);
    ylim(axes3,[3 6]);
    yticks([3 4 5 6]);
elseif number==82
    xlim(axes3,[-8 8]);
    xticks([-5 0 5]);
    ylim(axes3,[-7 1]);
    yticks([-7 -3 1]);
elseif number==83
    xlim(axes3,[-8 12]);
    xticks([-6 0 6 12]);
    ylim(axes3,[-0.03 0.005]);
    yticks([-0.02 0]);
elseif number==91
    xlim(axes3,[-10 10]);
    xticks([-10 -5 0 5 10]);
    ylim(axes3,[0 160]);
    yticks([0 80 160]);
elseif number==92
    xlim(axes3,[-10 10]);
    xticks([-10 -5 0 5 10]);
    ylim(axes3,[-80 80]);
    yticks([-80 0 80]); 
end
     
% Create 4 plots
plot(x,variationfreeenergy(:,1),'DisplayName',times(1,:),'LineWidth',3,'LineStyle',':',...
    'Color',[0 0.5647 0.6196]);

hold all
plot(x,variationfreeenergy(:,2),'DisplayName',times(2,:),'LineWidth',3,'LineStyle','--',...
    'Color',[0.4980 0.6902 0.0196]);

hold all
plot(x,variationfreeenergy(:,3),'DisplayName',times(3,:),'LineWidth',3,'LineStyle','-.',...
    'Color',[0.87058824300766 0.490196079015732 0]);

hold all
plot(x,variationfreeenergy(:,4),'DisplayName',times(4,:),'LineWidth',3,'LineStyle','-',...
    'Color',[0.4863 0.0784 0.3020]);

% Set the remaining axes properties
box(axes3,'on');
set(axes3,'FontSize',16,'TickLabelInterpreter','latex','XGrid','on','YGrid',...
    'on');

% Create xlabel
xlabel({'$$x$$'},'FontSize',23,'Interpreter','latex');

% Create legend
if number==82
    legend3 = legend(axes3,'show','Location','southeast');
else
    legend3 = legend(axes3,'show','Location','northeast');
end
set(legend3,'Interpreter','latex','FontSize',21);

% Create textarrow for the ylabel
annotation(figure3,'textarrow',[0.033 0.1],...
    [0.584 0.865459627329192],'TextEdgeColor','none',...
    'FontSize',23,...
    'String',{'$$\frac{\delta \mathcal{F}}{\delta \rho}$$'},...
    'HeadStyle','none',...
    'LineStyle','none','Interpreter','latex');

% Save the figure
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
saveas(gcf, [pwd sprintf('/figures/varfreeenergy-%d',number)], 'pdf')%Save figure

%--------------------------------------------------------------------------
% Forth plot: temporal evolution of the total and free energy
%-------------------------------------------------------------------------- 

% Create figure
figure4 = figure('PaperOrientation','landscape','pos',[10 10 560 450],'PaperSize',[20 13]);

% Create axes
axes4 = axes('Parent',figure4);
hold(axes4,'on');

% Customize axes for each example
if number==1
    xlim(axes4,[0 7]);
    xticks([0 5]);
    ylim(axes4,[-2  0]);
    yticks([-2 -1 0]);
elseif number==2
    xlim(axes4,[0 7]);
    xticks([0 5]);
    ylim(axes4,[-2  0]);
    yticks([-2 -1 0]);
elseif number==3
    xlim(axes4,[0 20]);
    xticks([0 10 20]);
    ylim(axes4,[-2  0]);
    yticks([-2 -1 0]);
elseif number==4
    xlim(axes4,[0 6]);
    xticks([0 3 6]);
    ylim(axes4,[0  4.9]);
    yticks([0 2 4]);
elseif number==5
    xlim(axes4,[0 3]);
    xticks([0 1.5 3]);
    ylim(axes4,[-1.925  -1.885]);
    yticks([-1.92 -1.89]);
elseif number==61 || number==62 || number==63
    xlim(axes4,[0 6]);
    xticks([0 4]);
    ylim(axes4,[-3  10]);
    yticks([0 5 10]);
elseif number==7
    xlim(axes4,[0 7]);
    xticks([0 5]);
    ylim(axes4,[-2  0]);
    yticks([-2 -1 0]);
elseif number==81
    xlim(axes4,[0 20]);
    xticks([0 15]);
    ylim(axes4,[1.95  2.25]);
    yticks([1.95 2.1  2.25]);
elseif number==82
    xlim(axes4,[0 200]);
    xticks([0 100 200]);
    ylim(axes4,[-1.5  1.5]);
    yticks([-1.5 0 1.5]);
elseif number==83
    xlim(axes4,[0 1200]);
    xticks([0 600 1200]);
    ylim(axes4,[-0.021  -0.013]);
    yticks([-0.021 -0.017 -0.013]);
elseif number==91
    xlim(axes4,[0 6]);
    xticks([0 3 6]);
    ylim(axes4,[55  90]);
    yticks([60 75 90]);
elseif number==92
    xlim(axes4,[0 16]);
    xticks([0 8 16]);
    ylim(axes4,[-15  15]);
    yticks([-15 0 15]);
end
   
% Create 4 plots
plot(t,totalenergy,'DisplayName','$$E^\Delta$$','LineWidth',3,...
    'Color',[0 0.5647 0.6196]);

hold all
plot(t,freeenergy,'DisplayName','$$F^\Delta$$','LineWidth',3,'LineStyle',':',...
    'Color',[0.87058824300766 0.490196079015732 0]);

% Set the remaining axes properties
box(axes4,'on');
set(axes4,'FontSize',16,'TickLabelInterpreter','latex','XGrid','on','YGrid',...
    'on');

% Create xlabel
xlabel({'$$t$$'},'FontSize',23,'Interpreter','latex');

% Create legend
legend4 = legend(axes4,'show','Location','northeast');
set(legend4,'Interpreter','latex','FontSize',21);

% Save the figure
set(gcf, 'PaperPosition', [0 0 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
saveas(gcf, [pwd sprintf('/figures/energies-%d',number)], 'pdf')%Save figure
end


