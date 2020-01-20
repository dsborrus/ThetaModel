data = O.spikes;

[PKS,LOCS,w] = findpeaks(data,O.t/1000,'minpeakwidth',.25,'minpeakheight',1);

gaps = diff(LOCS);

figure; hold on;
title('preintervVSburstsize')
plot(gaps,PKS(2:end),'ro');
xlabel('preceding interval'); ylabel('burst amplitude');

% 
% figure; hold on;
% title('preintervVSburstsize')
% plot(gaps,w(2:end),'ro');
% xlabel('preceding interval'); ylabel('burst area?');