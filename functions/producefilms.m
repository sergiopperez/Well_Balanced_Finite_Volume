function producefilms(n,x,U,t,numberexample)


writerObj = VideoWriter([pwd sprintf('/films/example-%d',numberexample)],'Motion JPEG AVI');
writerObj.FrameRate = 16;

writerObj.Quality = 100;
open(writerObj);

timescreenshot=0;

% Define for each frame
if numberexample==1 
    timebetweenframes=0.09;
    tmaxfilm=13;
elseif numberexample==2 
    timebetweenframes=0.09;
    tmaxfilm=13;
elseif numberexample==3 
    timebetweenframes=0.25;
    tmaxfilm=70;
elseif numberexample==4 
    timebetweenframes=0.075;
    tmaxfilm=10;
elseif numberexample==5 
    timebetweenframes=0.04;
    tmaxfilm=3;
elseif numberexample==61 || numberexample==62 || numberexample==63 
    timebetweenframes=0.07;
    tmaxfilm=10;
elseif numberexample==7 
    timebetweenframes=0.09;
    tmaxfilm=13;
elseif numberexample==81 
    timebetweenframes=0.16;
    tmaxfilm=25;
elseif numberexample==82 
    timebetweenframes=0.35;
    tmaxfilm=65;
elseif numberexample==83 
    timebetweenframes=2.5;
    tmaxfilm=500;
elseif numberexample==91 
    timebetweenframes=0.04;
    tmaxfilm=6;
elseif numberexample==92 
    timebetweenframes=0.08;
    tmaxfilm=15;
end
numbertime=1;
for i=1:length(t)
    
    
    if t(i)>=timescreenshot
        
        
        
        %q=str2num(num2str(T(1,i),'%10.2e'));
        %q=str2num(num2str(t(i,1),'%10.3f'));
        %q=num2str(t(i,1),'%10.4f');
        q=num2str(t(i,1),'t=%.2f');
        %q=num2str(t(i,1),'%05.5g');
%         if i==1
%             q='0.0000'
%         end
            
        %figure('Visible','Off')
        
        writerObj=createfigure_density_momentum(x', U(1:n,i), U(n+1:2*n,i),q,numberexample,numbertime,writerObj);
        
%         ax = gca;
% ax.Units = 'pixels';
% pos = ax.Position;
% ti = ax.TightInset;
% rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
% frame= getframe(ax,rect);
%frame= getframe(figure1);

       % frame = getframe(gcf);%, [0 0 849 720])%436 343]);
        
%        writeVideo(writerObj,frame);
        
        timescreenshot=timescreenshot+timebetweenframes;

        numbertime=numbertime+1;
        
        if t(i)>tmaxfilm
        break
        end
        
    end
end

close(writerObj);

end

