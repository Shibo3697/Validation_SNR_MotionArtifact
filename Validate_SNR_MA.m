clear; clc; close all;

%% path
dataFolder = 'C:\Users\admin\Desktop\DOTBallsqueezing'; 
filePattern = fullfile(dataFolder, '**', '*.snirf');
snirfFiles = dir(filePattern);


%% Data Containers
results = struct();
results.Rest.SNR = []; 
results.Rest.MotionRatio = [];
results.Task.SNR = []; 
results.Task.MotionRatio = [];
results.Motion.SNR = []; 
results.Motion.MotionRatio = [];

%% Main Processing Loop
h = waitbar(0, 'Initializing analysis...');

for k = 1 : length(snirfFiles)
    % Update Progress Bar
    waitbar(k/length(snirfFiles), h, sprintf('Processing file %d of %d...', k, length(snirfFiles)));
    
    fullFileName = fullfile(snirfFiles(k).folder, snirfFiles(k).name);
    snirfObj = SnirfClass(fullFileName);

    
    % ---  SNR ---
    d = snirfObj.data.dataTimeSeries;  
    channel_SNR = 20*log10(mean(d)./std(d));
    avg_SNR = mean(channel_SNR, 'omitnan'); 
    
    % --- Motion Artifact Ratio ---
    % [tInc, tIncCh] = hmrR_MotionArtifactByChannel(hmrR_Intensity2OD(snirfObj.data), snirfObj.probe, [], [], [], 0.5, 1.0, 10, 2);
    % motion_ratio_per_channel = (size(tIncCh{1}, 1) - sum(tIncCh{1}, 1)) ./ size(tIncCh{1}, 1);
    % avg_MR = mean(motion_ratio_per_channel, 'omitnan');
    DOD = hmrR_Intensity2OD(snirfObj.data);
    [tInc, tIncCh]=hmrR_MotionArtifactByChannel(DOD,snirfObj.probe, [], [], [], 0.5, 0.5, 20, 5);
    tIncCh_bin = tIncCh{1} > 0;
    motion_ratio_per_channel = 1 - mean(tIncCh_bin, 1);
    avg_MR = mean(motion_ratio_per_channel, 'omitnan');
    
    % --- Categorize ---
    if contains(snirfFiles(k).name, 'task-Rest', 'IgnoreCase', true)
        results.Rest.SNR(end+1) = avg_SNR;
        results.Rest.MotionRatio(end+1) = avg_MR;
    elseif contains(snirfFiles(k).name, 'task-BallSqueezing', 'IgnoreCase', true)
        results.Task.SNR(end+1) = avg_SNR;
        results.Task.MotionRatio(end+1) = avg_MR;
    elseif contains(snirfFiles(k).name, 'task-Motion', 'IgnoreCase', true)
        results.Motion.SNR(end+1) = avg_SNR;
        results.Motion.MotionRatio(end+1) = avg_MR;
    end
end
close(h); % Close progress bar

%% Calculate Statistics
% SNR
mean_Rest_SNR = mean(results.Rest.SNR); std_Rest_SNR = std(results.Rest.SNR);
mean_Task_SNR = mean(results.Task.SNR); std_Task_SNR = std(results.Task.SNR);
mean_Motion_SNR = mean(results.Motion.SNR); std_Motion_SNR = std(results.Motion.SNR);

% Motion Ratio (Convert to %)
mean_Rest_MR = mean(results.Rest.MotionRatio)*100; std_Rest_MR = std(results.Rest.MotionRatio)*100;
mean_Task_MR = mean(results.Task.MotionRatio)*100; std_Task_MR = std(results.Task.MotionRatio)*100;
mean_Motion_MR = mean(results.Motion.MotionRatio)*100; std_Motion_MR = std(results.Motion.MotionRatio)*100;

%% Print Result
fprintf('\n======================================================\n');
fprintf('       TECHNICAL VALIDATION RESULTS SUMMARY           \n');
fprintf('======================================================\n');

fprintf('1. RESTING STATE (ses-01)\n');
fprintf('   Mean SNR: %.2f +/- %.2f dB\n', mean_Rest_SNR, std_Rest_SNR);
fprintf('   Artifact Ratio: %.2f%% +/- %.2f%%\n', mean_Rest_MR, std_Rest_MR);

fprintf('\n2. BALL SQUEEZING TASK (ses-02)\n');
fprintf('   Mean SNR: %.2f +/- %.2f dB\n', mean_Task_SNR, std_Task_SNR);
fprintf('   Artifact Ratio: %.2f%% +/- %.2f%%\n', mean_Task_MR, std_Task_MR);

fprintf('\n3. INSTRUCTED MOTION (ses-03)\n');
fprintf('   Mean SNR: %.2f +/- %.2f dB\n', mean_Motion_SNR, std_Motion_SNR);
fprintf('   Artifact Ratio: %.2f%% +/- %.2f%%\n', mean_Motion_MR, std_Motion_MR);
fprintf('------------------------------------------------------\n');

%% Visualization
figure('Name', 'Technical Validation', 'Color', 'w', ...
       'Position', [100, 100, 1200, 450]);

labels = {'Rest', 'Ball-Squeeze', 'Motion'};

%% Plot A
subplot(1,2,1)

data_SNR = [results.Rest.SNR(:); ...
            results.Task.SNR(:); ...
            results.Motion.SNR(:)];
        
group_SNR = [ones(numel(results.Rest.SNR),1); ...
             2*ones(numel(results.Task.SNR),1); ...
             3*ones(numel(results.Motion.SNR),1)];

boxplot(data_SNR, group_SNR, 'Labels', labels);
title('A. Signal Quality', 'FontSize', 14);
ylabel('Signal-to-Noise Ratio (dB)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'FontWeight', 'bold');
set(findobj(gca,'type','line'), 'LineWidth', 1.5);
grid on;


%% Plot B
subplot(1,2,2)

data_MR = 100 * [results.Rest.MotionRatio(:); ...
                 results.Task.MotionRatio(:); ...
                 results.Motion.MotionRatio(:)];

group_MR = [ones(numel(results.Rest.MotionRatio),1); ...
            2*ones(numel(results.Task.MotionRatio),1); ...
            3*ones(numel(results.Motion.MotionRatio),1)];

boxplot(data_MR, group_MR, 'Labels', labels);
title('B. Motion Artifact', 'FontSize', 14);
ylabel('Motion Artifact Ratio (%)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'FontWeight', 'bold');
set(findobj(gca,'type','line'), 'LineWidth', 1.5);
grid on;
ylim([-0.5, max(data_MR)*1.2]);

% t-test
%  Motion vs resting
[h1, p_val_rest] = ttest2(results.Motion.MotionRatio, results.Rest.MotionRatio, 'Vartype', 'unequal');
%  Motion vs Task 
[h2, p_val_task] = ttest2(results.Motion.MotionRatio, results.Task.MotionRatio, 'Vartype', 'unequal');

fprintf('Motion vs Rest: p = %.5f\n', p_val_rest);
fprintf('Motion vs Task: p = %.5f\n', p_val_task);