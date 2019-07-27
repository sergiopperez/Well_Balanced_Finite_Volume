function createfigure3D(x,U,t,number,index)
n=length(U(:,1))/2;

% Create figure
figure1 = figure('PaperOrientation','landscape','pos',[10 10 560 450],'PaperSize',[20 13]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% select number of lines
nlines=25;
tstep=t(index)/(nlines-1);
tnextline=0;
% for i=1:length(t)
%     if t(i)>=tnextline && t(i)<=t(index)
%      plot3(x,t(i)*ones(n,1),U(1:n,i)+0.002,'LineWidth',3,'LineStyle','-',...
%     'Color',[0 0.5647 0.6196]);
% tnextline=tnextline+tstep;
%     end
% end

% Collect only some times for the surface plot


t2=t([linspace(1,5001,51)']);




surface=surf(repmat(x,1,length(t2)),repmat(t2,1,length(x))',U(1:n,[linspace(1,5001,51)']));
set(surface,'FaceColor',[1 1 1],'MeshStyle','column','edgecolor',[0 0.5647 0.6196],'LineWidth',1,'FaceAlpha',1);

%shading interp

% mycolors = [1 1 1];
% colormap(mycolors);
view([20 60]);
if number==83
    xlim(axes1,[-8 12]);
    xticks([-6 0 6 12]);
    ylim(axes1,[0 400]);
    yticks([0 200 400]);
    zlim(axes1,[0 0.25]);
    zticks([0 0.2]);   
end

% Create xlabel

set(axes1,'FontSize',16,'TickLabelInterpreter','latex');

annotation(figure1,'textarrow',[0.033 0.1],...
    [0.43 0.865459627329192],'TextEdgeColor','none',...
    'FontSize',23,...
    'FontName','TeX Gyre Bonum',...
    'String',{'$$\rho$$'},...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'LineStyle','none');

annotation(figure1,'textarrow',[0.373 0.1],...
    [0.1 0.865459627329192],'TextEdgeColor','none',...
    'FontSize',23,...
    'String',{'$$x$$'},...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'LineStyle','none');

annotation(figure1,'textarrow',[0.923 0.1],...
    [0.30 0.895459627329192],'TextEdgeColor','none',...
    'FontSize',23,...
    'String',{'$$t$$'},...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'LineStyle','none');

set(gcf, 'PaperPosition', [0 1 15 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 15]); %Set the paper to have width 5 and height 5.
saveas(gcf, sprintf('density3D-%d',number), 'pdf') %Save figure
end