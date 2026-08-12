clear; clc; close all;

dataRootDir = '/Users/zhoushibo/Desktop/datapaper_data/sub-182';
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
% extract information
motionIdx = find(contains({fileList.name}, 'task-Motion'), 1);
motionPath = fullfile(fileList(motionIdx).folder, fileList(motionIdx).name);
snirf = SnirfClass(motionPath);

% Marker name
markerLabels = {'l-tilt', 'r-tilt', 'hello', 'brows', 'd-tilt', 'u-tilt'};
markerColors = [0.8 0.2 0.2;   % 1 - tilt head left
                0.2 0.2 0.8;   % 2 - tilt head right
                0.2 0.7 0.2;   % 3 - say hello
                0.9 0.5 0.0;   % 4 - lifting eyebrows
                0.5 0.0 0.5;   % 5 - tilt head down
                0.0 0.6 0.6];  % 6 - tilt head up

% extract event 
evTimes = [];
evTypes = [];
for i = 1:length(snirf.stim)
    markerID = str2double(snirf.stim(i).name);  % name not index
    onsets = snirf.stim(i).data(:,1);
    mask = onsets <= 200;
    evTimes = [evTimes; onsets(mask)];
    evTypes = [evTypes; repmat(markerID, sum(mask), 1)];
end
[evTimes, sortIdx] = sort(evTimes);
evTypes = evTypes(sortIdx);

%% Plot - Motion Raw Intensity with labeled markers
figure('Color', 'w', 'Position', [150, 50, 1000, 400]);
hold on;

for j = 1:length(evTimes)
    patch([evTimes(j), evTimes(j)+1, evTimes(j)+1, evTimes(j)], ...
        [0.55, 0.55, 0.7, 0.7], [0.9 0.9 0.9], ...
        'EdgeColor', 'none', 'HandleVisibility', 'off');
end

h_data = plot(snirf.data.time, snirf.data.dataTimeSeries(:, channelIdx), ...
    'Color', '#0072BD', 'LineWidth', 1.2);

for j = 1:length(evTimes)
    xline(evTimes(j), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    text(evTimes(j), 0.695, markerLabels{evTypes(j)}, ...
        'FontSize', 8, 'FontWeight', 'bold', ...
        'Color', markerColors(evTypes(j),:), ...
        'HorizontalAlignment', 'center', ...
        'Rotation', 45);
end

xlim([0, 200]); ylim([0.55 0.7]);
title('Motion Creation Task', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Intensity', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 12);
legend(h_data, 'Intensity', 'Location', 'northeast');

%% Accelerometer with labeled markers
figure('Color', 'w', 'Position', [150, 50, 1000, 800]);

for s = 1:2
    subplot(2, 1, s);
    hold on;

    if s == 1
        idx = [1, 2, 3]; t_str = 'accelerometer\_Cz';
    else
        idx = [13, 14, 15]; t_str = 'accelerometer\_Pz';
    end

    for j = 1:length(evTimes)
        patch([evTimes(j), evTimes(j)+1, evTimes(j)+1, evTimes(j)], ...
            [-12, -12, 14, 14], [0.9 0.9 0.9], ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    % Accelerometer data
    t_axis = snirf.aux(idx(1)).time;
    h1 = plot(t_axis, snirf.aux(idx(1)).dataTimeSeries, 'r');
    h2 = plot(t_axis, snirf.aux(idx(2)).dataTimeSeries, 'g');
    h3 = plot(t_axis, snirf.aux(idx(3)).dataTimeSeries, 'b');

    % marker
    for j = 1:length(evTimes)
        xline(evTimes(j), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        text(evTimes(j), 11, markerLabels{evTypes(j)}, ...
            'FontSize', 8, 'FontWeight', 'bold', ...
            'Color', markerColors(evTypes(j),:), ...
            'HorizontalAlignment', 'center', ...
            'Rotation', 45);
    end

    title(t_str, 'FontSize', 16, 'FontWeight', 'bold');
    if s == 2
        xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
    end
    xlim([0, 200]); ylim([-12 14]); grid on; box off;
    legend([h1, h2, h3], 'X', 'Y', 'Z', 'Location', 'northeast');
end

axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.08, 0.5, 'Acceleration (m/s^2)', 'FontSize', 16, 'FontWeight', 'bold', ...
    'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
