% time test

repeat = 10;

% OUTERt = tic;
% for k = 1:repeat
%     synctheta_v3;
%     k
% end
% tOUTER1 = toc(OUTERt);

tic
OUTERt = tic;
for k = 1:repeat
    synctheta_v3B; clearvars -except OUTERt k repeat
    k
end
tOUTER2 = toc(OUTERt);

OUTERt = tic;
for k = 1:repeat
    synctheta_v3C; clearvars -except OUTERt k tOUTER2 repeat
    k
end
tOUTER3 = toc(OUTERt);
disp(['Average time for case 1 was ' mat2str(tOUTER2/repeat)])
disp(['Average time for case 2 waw ' mat2str(tOUTER3/repeat)])