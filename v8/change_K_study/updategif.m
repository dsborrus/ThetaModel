function [P,P2,P3] = updategif(curtheta,pmin,pmax,Isummed,gf,gffilename,j,P,P2,P3,raster,spikes,t,vin,istate)
% meant to be used in synctheta  

D = vin(1);
isig1 = vin(2);
iu1 = vin(3);
isig2 = vin(4);
iu2 = vin(5);
N = vin(6);
prob = vin(7);
tauavg = vin(8);
mgain = vin(9);
taum = vin(10);
nrise = vin(11);
taun = vin(12);
sigain = vin(13);
tausi = vin(14);
tautheta=vin(15);

pres = 500;
parray = linspace(-.5,.5,pres);
theIs = NaN(N,pres);
theIs = (1/tautheta)*(1-cos(parray))+Isummed'*(1+cos(parray));

subplot(6,6,[1 21]); hold on;
delete(P); delete(P2); delete(P3);
P2 = plot(parray,theIs,'r','linewidth',0.1); hold on;
P = plot(parray,mean(theIs),'b','linewidth',2);
P3 = plot(curtheta,0,'ko','markersize',12);

tw = j-20*2000:j;

% network %
subplot(6,6,[4 6]); hold on;
cla;
xlim([t(tw(1))/1000 t(tw(end))/1000])
plot(t(tw)/1000,(spikes(tw)./tauavg)*1000)

% raster %
subplot(6,6,[10 12]); hold on; cla;
xlim([t(tw(1))/1000 t(tw(end))/1000])
for i = 1:N
    plot(t(tw)/1000,raster((tw),i)*i,'.k'); hold on;
end


% %
drawnow

frame = getframe(gf);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

if j==60001
    imwrite(imind,cm,gffilename,'gif','Loopcount',inf);
else
    imwrite(imind,cm,gffilename,'gif','WriteMode','append');
end
