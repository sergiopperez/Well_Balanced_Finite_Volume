%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function writerObj=createfigure_density_momentum(X1, Y1, Y2,q,numberexample,numbertime,writerObj)


%--------------------------------------------------------------------------
% First plot: density profile
%-------------------------------------------------------------------------- 

figure1 = figure('Color',...
    [0.803921580314636 0.878431379795074 0.968627452850342],'pos',[0 0 1100 1000], 'visible', 'off','PaperOrientation','landscape');
% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.572446555819477 0.775 0.311163895486936],...
    'FontName','TeX Gyre Bonum','FontSize',20,'TickLabelInterpreter', 'latex');


% Customize axes for each example
if numberexample==1
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.6]);
    yticks([0 0.3 0.6]);
elseif numberexample==2
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.9]);
    yticks([0 0.3 0.6 0.9]);
elseif numberexample==3
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 1.2]);
    yticks([0 0.5 1.0]);
elseif numberexample==4
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 1.0]);
    yticks([0 0.5 1.0]);
elseif numberexample==5
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.5]);
    yticks([0 0.2 0.4]); 
elseif numberexample==61 
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.8]);
    yticks([0 0.4 0.8]);
elseif numberexample==62 
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.9]);
    yticks([0 0.45 0.9]);
elseif numberexample==63 
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 1.0]);
    yticks([0 0.5 1.0]);
elseif numberexample==7
    xlim(axes1,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.7]);
    yticks([0 0.3 0.6]);
elseif numberexample==81
    xlim(axes1,[-8 8]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.4]);
    yticks([0 0.2 0.4]);
elseif numberexample==82
    xlim(axes1,[-8 8]);
    xticks([-5 0 5]);
    ylim(axes1,[0 0.4]);
    yticks([0 0.2 0.4]);
elseif numberexample==83
    xlim(axes1,[-8 12]);
    xticks([-6 0 6 12]);
    ylim(axes1,[0 0.24]);
    yticks([0 0.1 0.2]);
elseif numberexample==91
    xlim(axes1,[-6 6]);
    xticks([-6 -3 0 3 6]);
    ylim(axes1,[0 2]);
    yticks([0 1 2]);
elseif numberexample==92
 xlim(axes1,[-6 6]);
    xticks([-6 -3 0 3 6]);
    ylim(axes1,[0 2]);
    yticks([0 1 2]);  
end

box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(X1,Y1,'Parent',axes1,'MarkerSize',12,'Marker','.','LineWidth',3,...
    'LineStyle','none');

% Create xlabel
xlabel({'$$x$$'},'FontSize',26,'FontName','TeX Gyre Bonum','interpreter','latex');

% Create title
title(['EVOLUTION OF THE DENSITY'],...
    'FontWeight','demi',...
    'FontSize',26,...
    'FontName','TeX Gyre Bonum','interpreter','latex');


%--------------------------------------------------------------------------
% Second plot: momentum profile
%-------------------------------------------------------------------------- 



% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.13 0.0983449883449883 0.775 0.311163895486936],...
    'FontName','TeX Gyre Bonum','FontSize',20,'TickLabelInterpreter', 'latex');

% Customize axes for each example
if numberexample==1
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.25 0.25]);
    yticks([-0.25 0 0.25]);
elseif numberexample==2
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.35 0.35]);
    yticks([-0.3 0 0.3]);
elseif numberexample==3
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.7 0.7]);
    yticks([-0.7 0 0.7]);
elseif numberexample==4
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.5 0.5]);
    yticks([-0.5 0 0.5]);
elseif numberexample==5
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
%     ylim(axes2,[0 0.09]);
%     yticks([0 0.08]);
    ylim(axes2,[0 0.4]);
    yticks([0 0.2 0.4]); 
elseif numberexample==61  
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.4 0.4]);
    yticks([-0.4 0 0.4]);
elseif numberexample==62 
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.5 0.5]);
    yticks([-0.5 0 0.5]);
elseif numberexample==63
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.9 0.9]);
    yticks([-0.9 0 0.9]);
elseif numberexample==7
    xlim(axes2,[-5 5]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.4 0.4]);
    yticks([-0.4 0 0.4]);
elseif numberexample==81
    xlim(axes2,[-8 8]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.04 0.04]);
    yticks([-0.04 0 0.04]);
elseif numberexample==82
    xlim(axes2,[-8 8]);
    xticks([-5 0 5]);
    ylim(axes2,[-0.08 0.08]);
    yticks([-0.08 0 0.08]);
elseif numberexample==83
    xlim(axes2,[-8 12]);
    xticks([-6 0 6 12]);
    ylim(axes2,[-0.005 0.005]);
    yticks([-0.005 0 0.005])
elseif numberexample==91
    xlim(axes2,[-6 6]);
    xticks([-6 -3 0 3 6]);
    ylim(axes2,[-1 1]);
    yticks([-1 0 1]);
elseif numberexample==92
 xlim(axes2,[-6 6]);
    xticks([-6 -3 0 3 6]);
    ylim(axes2,[-1.3 1.3]);
    yticks([-1 0 1]);
end

box(axes2,'on');
hold(axes2,'all');

% Create plot
if numberexample~=5
plot(X1,Y2,'Parent',axes2,'MarkerSize',10,'Marker','.','LineWidth',2,...
    'LineStyle','none');
elseif numberexample==5
    plot(X1,Y2./Y1,'Parent',axes2,'MarkerSize',10,'Marker','.','LineWidth',2,...
    'LineStyle','none');
end
% Create title
if numberexample~=5
title({'EVOLUTION OF THE MOMENTUM'},'FontWeight','demi','FontSize',26,...
    'FontName','TeX Gyre Bonum','interpreter','latex');
elseif numberexample==5
    title({'EVOLUTION OF THE VELOCITY'},'FontWeight','demi','FontSize',26,...
    'FontName','TeX Gyre Bonum','interpreter','latex');
end
% Create xlabel
xlabel('$$x$$','FontSize',26,'FontName','TeX Gyre Bonum','interpreter','latex');

% Create ylabel
ylabel({''},'FontName','TeX Gyre Bonum');

% Create textbox
annotation(figure1,'textbox',...
    [0.05 0.455631823688853 0.14140639336261594 0.05124977150507],...
    'String',q,...
    'HorizontalAlignment','center','FontSize',21,...
    'VerticalAlignment','middle',...
    'FontName','TeX Gyre Bonum','interpreter','latex',...
    'EdgeColor',[0 0.447058826684952 0.74117648601532],...
    'LineWidth',1,...
    'BackgroundColor',[0.87058824300766 0.921568632125854 0.980392158031464]);

% Create textarrow
annotation(figure1,'textarrow',[0.068 0.166649810366625],...
    [0.75 0.865459627329192],'TextEdgeColor','none',...
    'FontSize',26,...
    'FontName','TeX Gyre Bonum',...
    'interpreter','latex',...
    'String',{'$$\rho$$'},...
    'HeadStyle','none',...
    'LineStyle','none');

% Create textarrow
if numberexample~=5
annotation(figure1,'textarrow',[0.07 0.16],...
    [0.27 0.396701863354036],'TextEdgeColor','none',...
    'FontName','TeX Gyre Bonum',...
    'interpreter','latex',...
    'FontSize',26,...
    'String',{'$$\rho u$$'},...
    'HeadStyle','none',...
    'LineStyle','none');
elseif numberexample==5
    annotation(figure1,'textarrow',[0.07 0.16],...
    [0.27 0.396701863354036],'TextEdgeColor','none',...
    'FontName','TeX Gyre Bonum',...
    'interpreter','latex',...
    'FontSize',26,...
    'String',{'$$u$$'},...
    'HeadStyle','none',...
    'LineStyle','none');
    
end
print([pwd sprintf('/screenshots/example-%d-frame-%g',numberexample,numbertime)],'-dpdf','-fillpage')

frame= getframe(figure1);
writeVideo(writerObj,frame);




