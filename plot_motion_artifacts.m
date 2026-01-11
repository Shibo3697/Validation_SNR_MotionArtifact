clc;
clear;

filename = 'sub-176_ses-03_task-Motion_nirs.snirf';
snirf = SnirfClass(filename);

% find a channel
[~, ch] = max(std(snirf.data.dataTimeSeries));

figure('Color', 'w');


for i = 1:length(snirf.stim)
    onset = snirf.stim(i).data(1, 1);
    
    % time window: 2s before to 8s after onset
    idx = find(snirf.data.time > (onset - 2) & snirf.data.time < (onset + 8));

    subplot(2, 3, i);
    plot(snirf.data.time(idx) - onset, snirf.data.dataTimeSeries(idx, ch), 'k'); 
    hold on;
    plot([0, 0], ylim, 'r--'); 
    title(snirf.stim(i).name, 'Interpreter', 'none'); 
    axis tight; 
end