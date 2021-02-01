
close all
clearvars;
% Load two example transects
load transect

% Separate bathymetric data  representing multiscale bedforms
% Example 1
NF = 15;
cutoff_slope = 0.03;
loess_br_fit1 = loess_br(transect.x1,transect.y1,NF,cutoff_slope);

figure; 
subplot(3,1,1); hold on; plot(transect.x1, transect.y1, 'k'); 
    plot(transect.x1, loess_br_fit1,'r');
    xlim([min(transect.x1) max(transect.x1)]); ylabel('y [m]'); xlabel('x [m]');
subplot(3,1,2); hold on; plot(transect.x1, transect.y1-loess_br_fit1,'k');
    plot(transect.x1,zeros(length(transect.x1),1),'r')
    xlim([min(transect.x1) max(transect.x1)]); ylabel('y [m]'); xlabel('x [m]');
subplot(3,1,3); hold on; plot(transect.x1, transect.y1-loess_br_fit1,'k');
    plot(transect.x1,zeros(length(transect.x1),1),'r')
xlim([420 480]); ylabel('y [m]'); xlabel('x [m]');

% Example 2
NF = 20;
cutoff_slope = 0.03;
loess_br_fit2 = loess_br(transect.x2,transect.y2,NF,cutoff_slope);

figure; 
subplot(3,1,1); hold on; plot(transect.x2, transect.y2, 'k'); 
    plot(transect.x2, loess_br_fit2,'r');
    xlim([min(transect.x2) max(transect.x2)]); ylabel('y [m]'); xlabel('x [m]');
subplot(3,1,2); hold on; plot(transect.x2, transect.y2-loess_br_fit2,'k');
    plot(transect.x2,zeros(length(transect.x2),1),'r')
    xlim([min(transect.x2) max(transect.x2)]); ylabel('y [m]'); xlabel('x [m]');
subplot(3,1,3); hold on; plot(transect.x2, transect.y2-loess_br_fit2,'k');
    plot(transect.x2,zeros(length(transect.x2),1),'r')
xlim([420 480]); ylabel('y [m]'); xlabel('x [m]');

