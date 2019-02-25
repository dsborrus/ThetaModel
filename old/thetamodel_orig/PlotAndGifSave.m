% Plot and save
% Does all plotting and gif saving

hold off
%plot(exp(i*phase),'k-'); hold on;
%plot(1.4*exp(i*phase),'k-'); hold on;
plot(exp(-i*pi/2)*exp(i*theta(j-nplot+2:j+1,1:n1)).*(ones(nplot,1)*(1:n1))/n1,'go:'); hold on;
plot(exp(-i*pi/2)*exp(i*theta(j-nplot+2:j+1,n1+1:end)).*(ones(nplot,1)*(1:n2))/n2,'ro:');
plot([0 0],[0 1.05],'b:');
%plot(rplot(j-nplot:j+1),'ro')

%title(['t = ' num2str(floor(dt*j)) ])

axis equal
axis([ -1 1 -1 1.3 ])
axis off
box off
plot(-1+2*t/t(end),1.1+spikes/40,'b');


%title([ 'n=' num2str(n1) '&' num2str(n2) ',p=' num2str(prob) ', i=' num2str(iu1) '&' num2str(iu2) ', D=' num2str(D) ', \tau=' num2str(tauavg)])

drawnow

if 1
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if firstframe == 1;
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0);
        firstframe = 0;
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
    end
end