clear; clc; close all;

dataRootDir = 'C:\Users\t243f765\Desktop\DOTBallsqueezing\sub-182';
fileList = dir(fullfile(dataRootDir, '**', '*.snirf'));
 
timeWindow = [0,100];
channelIdx = 101;          % channel

% plot
figure('Color', 'w', 'Position', [150, 50, 1000, 800]); 

for i = 1:5
    filePath = fullfile(fileList(i).folder, fileList(i).name);
    [~, fileName, ~] = fileparts(fileList(i).name);
    snirfObj = SnirfClass(filePath); 
    d = snirfObj.data.dataTimeSeries(:, channelIdx); 
    t = snirfObj.data.time;
    subplot(5, 1, i);
    plot(t, d, 'Color', '#0072BD', 'LineWidth', 1.2); 
    
    xlim(timeWindow); 
    ylim([0.6 0.7]);
    
    title(fileName, 'FontSize', 12, 'Interpreter', 'none', 'FontWeight', 'bold');
    grid on;
    
    set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 12);
    if i == 3
        ylabel('Intensity', 'FontSize', 16, 'FontWeight', 'bold');
    end
    if i < 5
        set(gca, 'XTickLabel', []);
    end
end

xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
%% motion artifical
motionIdx = find(contains({fileList.name}, 'task-Motion'), 1);
motionPath = fullfile(fileList(motionIdx).folder, fileList(motionIdx).name);
snirf = SnirfClass(motionPath);

% Get all stimulus timing
evs = [];
for i = 1:length(snirf.stim)
    evs = [evs; snirf.stim(i).data(:,1)]; 
end
evs = evs(evs <= 200); 

% Plot
figure('Color', 'w', 'Position', [150, 50, 1000, 400]); 
hold on;
h_data = plot(snirf.data.time, snirf.data.dataTimeSeries(:, channelIdx), ...
              'Color', '#0072BD', 'LineWidth', 1.2);

for stim_t = evs' 
    p = patch([stim_t, stim_t+1, stim_t+1, stim_t], [0.55, 0.55, 0.7, 0.7], [0.9 0.9 0.9], ...
              'EdgeColor', 'none', 'HandleVisibility', 'off');
    shadow1 = xline(evs, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
end

xlim([0, 200]); 
title('Motion Creation Task', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Intensity', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 12);

legend([h_data, shadow1(1)], 'Intensity', 'Stimulus');
%% Accelerometer
figure('Color', 'w', 'Position', [150, 50, 1000, 800]); 

% plot
for s = 1:2
    subplot(2, 1, s); 
    hold on;
    for t = evs'
        patch([t, t+1, t+1, t], [-15, -15, 15, 15], [0.9 0.9 0.9], ...
              'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    shadow = xline(evs, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    
    if s == 1
        idx = [1, 2, 3]; t_str = 'accelerometer_Cz'; y_str = 'Acc 1';
    else
        idx = [13, 14, 15]; t_str = 'accelerometer_Pz'; y_str = 'Acc 3';
    end
    
    % Accelerometer data
    t_axis = snirf.aux(idx(1)).time;
    h1 = plot(t_axis, snirf.aux(idx(1)).dataTimeSeries, 'r'); 
    h2 = plot(t_axis, snirf.aux(idx(2)).dataTimeSeries, 'g'); 
    h3 = plot(t_axis, snirf.aux(idx(3)).dataTimeSeries, 'b');
    

    title(t_str, 'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'bold');
    if s == 2
        xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
    end
    xlim([0, 200]); ylim([-12 12]); grid on; box off;
    legend([h1, h2, h3, shadow(1)], 'X', 'Y', 'Z', 'Stimulus', 'Location', 'northeast');
end
axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.08, 0.5, 'Acceleration', 'FontSize', 16, 'FontWeight', 'bold', ...
     'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
