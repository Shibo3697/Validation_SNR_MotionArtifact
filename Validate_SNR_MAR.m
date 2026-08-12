clear; clc; close all;

%% function for plot
function addSigLine(x1, x2, y, h, label)
    line([x1 x2], [y y], 'LineWidth', 1.5, 'Color', 'k');
    line([x1 x1], [y-h y], 'LineWidth', 1.5, 'Color', 'k');
    line([x2 x2], [y-h y], 'LineWidth', 1.5, 'Color', 'k');
    text(mean([x1 x2]), y + h/4, label, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold');
end

%% Path
dataFolder   = 'C:\Users\t243f765\Desktop\Light Weight Multi-Distance fNIRS Dataset_';
filePattern  = fullfile(dataFolder, '**', '*.snirf');
snirfFiles   = dir(filePattern);

% Data Containers
results = struct();
results.Rest.SNR = [];   results.Rest.MR = [];
results.Task.SNR = [];   results.Task.MR = [];
results.Motion.SNR = []; results.Motion.MR = [];

%% Main Processing Loop
for k = 1 : numel(snirfFiles)
    fullFileName = fullfile(snirfFiles(k).folder, snirfFiles(k).name);
    snirfObj = SnirfClass(fullFileName);
    
    % SNR 
    d = snirfObj.data.dataTimeSeries;
    channel_SNR = 20*log10(mean(d)./std(d));
    
    % Motion Artifact Ratio 
    DOD = hmrR_Intensity2OD(snirfObj.data);
    [~, tIncCh] = hmrR_MotionArtifactByChannel(DOD, snirfObj.probe, [], [], [], 0.5, 0.5, 20, 5);
    tIncCh_bin = tIncCh{1} > 0;
    channel_MR = (1 - mean(tIncCh_bin, 1)) * 100; % MAR%
    
    % Categorize
fname = snirfFiles(k).name;
    if contains(fname, 'task-Rest', 'IgnoreCase', true)
        results.Rest.SNR(end+1, :) = channel_SNR;
        results.Rest.MR(end+1, :)  = channel_MR;
    elseif contains(fname, 'task-BallSqueezing', 'IgnoreCase', true)
        results.Task.SNR(end+1, :) = channel_SNR;
        results.Task.MR(end+1, :)  = channel_MR;
    elseif contains(fname, 'task-Motion', 'IgnoreCase', true)
        results.Motion.SNR(end+1, :) = channel_SNR;
        results.Motion.MR(end+1, :)  = channel_MR;
    end
end

% Average 3 blocks within subject for Ball-squeezing
results.Task.SNR = squeeze(mean(reshape(results.Task.SNR, 3, [], 200), 1));
results.Task.MR  = squeeze(mean(reshape(results.Task.MR,  3, [], 200), 1));

% Channel Distance Calculation
sdfile = load('probe_bilateral.SD', '-mat');
srcPos = sdfile.SD.SrcPos;
detPos = sdfile.SD.DetPos;
ml     = snirfObj.data(1).measurementList;

numChannels = length(ml);
dists = zeros(1, numChannels); % channel disatances

for i = 1:numChannels
    s = ml(i).sourceIndex;
    d = ml(i).detectorIndex;    
    dists(i) = norm(srcPos(s,:) - detPos(d,:)); % calculate distance, euclidean metric
end


threshold = 22;

idx_short = dists <= threshold; 
idx_long  = dists > threshold;

fprintf('shorter channel: %d , longer channel: %d\n', sum(idx_short), sum(idx_long));

% Split Data into Short/Long Groups for Analysis

% SNR Analysis
% Resting State
final_stats.Rest.SNR_Short = mean(results.Rest.SNR(:, idx_short), 2);
final_stats.Rest.SNR_Long  = mean(results.Rest.SNR(:, idx_long),  2);

% Ball Squeezing (Task)
final_stats.Task.SNR_Short = mean(results.Task.SNR(:, idx_short), 2);
final_stats.Task.SNR_Long  = mean(results.Task.SNR(:, idx_long),  2);

% Motion Creation
final_stats.Motion.SNR_Short = mean(results.Motion.SNR(:, idx_short), 2);
final_stats.Motion.SNR_Long  = mean(results.Motion.SNR(:, idx_long),  2);

% MAR Analysis
% Resting State
final_stats.Rest.MR_Short = mean(results.Rest.MR(:, idx_short), 2);
final_stats.Rest.MR_Long  = mean(results.Rest.MR(:, idx_long),  2);

% Ball Squeezing (Task)
final_stats.Task.MR_Short = mean(results.Task.MR(:, idx_short), 2);
final_stats.Task.MR_Long  = mean(results.Task.MR(:, idx_long),  2);

% Motion Creation
final_stats.Motion.MR_Short = mean(results.Motion.MR(:, idx_short), 2);
final_stats.Motion.MR_Long  = mean(results.Motion.MR(:, idx_long),  2);


condNames = {'Rest', 'Task', 'Motion'};
WithinStructure = table([1 2 3]', 'VariableNames', {'Conditions'});

%% ANOVA
% Short Channels (SNR)
fprintf('\n======================================================\n');
fprintf('[Analysis 1: SNR - Shorter Channels]\n');

% Normality Test
fprintf('--- Normality Test (Lilliefors) ---\n');

[h1,p1] = lillietest(final_stats.Rest.SNR_Short);   fprintf('  Rest:   p=%.3f (h=%d)\n', p1, h1);
[h2,p2] = lillietest(final_stats.Task.SNR_Short);   fprintf('  Task:   p=%.3f (h=%d)\n', p2, h2);
[h3,p3] = lillietest(final_stats.Motion.SNR_Short); fprintf('  Motion: p=%.3f (h=%d)\n', p3, h3);

t_snr_short = table(final_stats.Rest.SNR_Short, ...
                    final_stats.Task.SNR_Short, ...
                    final_stats.Motion.SNR_Short, ...
                    'VariableNames', condNames);
rm_snr_short = fitrm(t_snr_short, 'Rest-Motion~1', 'WithinDesign', WithinStructure);

% Mauchly Sphericity Test
mauchly_tbl = mauchly(rm_snr_short);
p_mauchly = mauchly_tbl.pValue;
fprintf('--- Mauchly Sphericity Test ---\n');
if p_mauchly < 0.05
    fprintf('  p = %.4f (<0.05). Violated. -> Look at [pValueGG] in ANOVA table.\n', p_mauchly);
else
    fprintf('  p = %.4f (>0.05). Met.      -> Look at [pValue] in ANOVA table.\n', p_mauchly);
end

tbl_snr_short = ranova(rm_snr_short);
fprintf('\n[SNR - Shorter Channels] ANOVA Result:\n');
disp(tbl_snr_short); 
fprintf('[SNR - Shorter] Pairwise Comparisons (Bonferroni):\n');
comp_snr_short = multcompare(rm_snr_short, 'Conditions', 'ComparisonType', 'bonferroni');
disp(comp_snr_short);

% Longer Channels (SNR)
fprintf('\n======================================================\n');
fprintf('[Analysis 2: SNR - Longer Channels]\n');

% Normality Test
fprintf('--- Normality Test (Lilliefors) ---\n');
[h1,p1] = lillietest(final_stats.Rest.SNR_Long);   fprintf('  Rest:   p=%.3f (h=%d)\n', p1, h1);
[h2,p2] = lillietest(final_stats.Task.SNR_Long);   fprintf('  Task:   p=%.3f (h=%d)\n', p2, h2);
[h3,p3] = lillietest(final_stats.Motion.SNR_Long); fprintf('  Motion: p=%.3f (h=%d)\n', p3, h3);

t_snr_long = table(final_stats.Rest.SNR_Long, ...
                   final_stats.Task.SNR_Long, ...
                   final_stats.Motion.SNR_Long, ...
                   'VariableNames', condNames);
rm_snr_long = fitrm(t_snr_long, 'Rest-Motion~1', 'WithinDesign', WithinStructure);

% Mauchly Sphericity Test
mauchly_tbl = mauchly(rm_snr_long);
p_mauchly = mauchly_tbl.pValue;
fprintf('--- Mauchly Sphericity Test ---\n');
if p_mauchly < 0.05
    fprintf('  p = %.4f (<0.05). Violated. -> Look at [pValueGG].\n', p_mauchly);
else
    fprintf('  p = %.4f (>0.05). Met.      -> Look at [pValue].\n', p_mauchly);
end

tbl_snr_long = ranova(rm_snr_long);
fprintf('[SNR - Longer Channels] ANOVA Result:\n');
disp(tbl_snr_long);
fprintf('[SNR - Longer] Pairwise Comparisons (Bonferroni):\n');
comp_snr_long = multcompare(rm_snr_long, 'Conditions', 'ComparisonType', 'bonferroni');
disp(comp_snr_long);

% Short Channels (MAR)
fprintf('\n======================================================\n');
fprintf('[Analysis 3: MAR - Shorter Channels]\n');

% Normality Test
fprintf('--- Normality Test (Lilliefors) ---\n');
[h1,p1] = lillietest(final_stats.Rest.MR_Short);   fprintf('  Rest:   p=%.3f (h=%d)\n', p1, h1);
[h2,p2] = lillietest(final_stats.Task.MR_Short);   fprintf('  Task:   p=%.3f (h=%d)\n', p2, h2);
[h3,p3] = lillietest(final_stats.Motion.MR_Short); fprintf('  Motion: p=%.3f (h=%d)\n', p3, h3);

t_mr_short = table(final_stats.Rest.MR_Short, ...
                   final_stats.Task.MR_Short, ...
                   final_stats.Motion.MR_Short, ...
                   'VariableNames', condNames);
rm_mr_short = fitrm(t_mr_short, 'Rest-Motion~1', 'WithinDesign', WithinStructure);

% Mauchly Sphericity Test
mauchly_tbl = mauchly(rm_mr_short);
p_mauchly = mauchly_tbl.pValue;
fprintf('--- Mauchly Sphericity Test ---\n');
if p_mauchly < 0.05
    fprintf('  p = %.4f (<0.05). Violated. -> Look at [pValueGG].\n', p_mauchly);
else
    fprintf('  p = %.4f (>0.05). Met.      -> Look at [pValue].\n', p_mauchly);
end

tbl_mr_short = ranova(rm_mr_short);
fprintf('[MAR - Shorter Channels] ANOVA Result:\n');
disp(tbl_mr_short);
fprintf('[MAR - Shorter] Pairwise Comparisons (Bonferroni):\n');
comp_mr_short = multcompare(rm_mr_short, 'Conditions', 'ComparisonType', 'bonferroni');
disp(comp_mr_short);

% Longer Channels (MAR)
fprintf('\n======================================================\n');
fprintf('[Analysis 4: MAR - Longer Channels]\n');

% Normality Test
fprintf('--- Normality Test (Lilliefors) ---\n');
[h1,p1] = lillietest(final_stats.Rest.MR_Long);   fprintf('  Rest:   p=%.3f (h=%d)\n', p1, h1);
[h2,p2] = lillietest(final_stats.Task.MR_Long);   fprintf('  Task:   p=%.3f (h=%d)\n', p2, h2);
[h3,p3] = lillietest(final_stats.Motion.MR_Long); fprintf('  Motion: p=%.3f (h=%d)\n', p3, h3);

t_mr_long = table(final_stats.Rest.MR_Long, ...
                  final_stats.Task.MR_Long, ...
                  final_stats.Motion.MR_Long, ...
                  'VariableNames', condNames);
rm_mr_long = fitrm(t_mr_long, 'Rest-Motion~1', 'WithinDesign', WithinStructure);

% Mauchly Sphericity Test
mauchly_tbl = mauchly(rm_mr_long);
p_mauchly = mauchly_tbl.pValue;
fprintf('--- Mauchly Sphericity Test ---\n');
if p_mauchly < 0.05
    fprintf('  p = %.4f (<0.05). Violated. -> Look at [pValueGG].\n', p_mauchly);
else
    fprintf('  p = %.4f (>0.05). Met.      -> Look at [pValue].\n', p_mauchly);
end

tbl_mr_long = ranova(rm_mr_long);
fprintf('[MAR - Longer Channels] ANOVA Result:\n');
disp(tbl_mr_long);
fprintf('[MAR - Longer] Pairwise Comparisons (Bonferroni):\n');
comp_mr_long = multcompare(rm_mr_long, 'Conditions', 'ComparisonType', 'bonferroni');
disp(comp_mr_long);

%% plot ANOVA
figure('Color', 'w', 'Position', [100, 100, 900, 700]);
labels = {'Resting', 'Ball-squeezing', 'Motion creation'};

% 1. SNR - Shorter Channels
subplot(2, 2, 1);
data = [final_stats.Rest.SNR_Short, final_stats.Task.SNR_Short, final_stats.Motion.SNR_Short];
boxplot(data, 'Labels', labels, 'Symbol', 'o');
title('SNR - Shorter Channels', 'FontSize', 16, 'FontWeight', 'bold'); 
ylabel('SNR (dB)', 'FontSize', 16, 'FontWeight', 'bold'); 
ylim([15 48]);
grid on;
% sig line
yMax = max(data(:));
y = yMax + 1;
h = 0.6;
addSigLine(1, 3, y, h, '**');   % Rest vs Motion
addSigLine(2, 3, y + 2.4, h, '***'); % ball vs Motion

% 2. SNR - Longer Channels
subplot(2, 2, 2);
data = [final_stats.Rest.SNR_Long, final_stats.Task.SNR_Long, final_stats.Motion.SNR_Long];
boxplot(data, 'Labels', labels, 'Symbol', 'o');
title('SNR - Longer Channels', 'FontSize', 16, 'FontWeight', 'bold'); 
ylabel('SNR (dB)', 'FontSize', 16, 'FontWeight', 'bold'); 
grid on;
ylim([15 48]);
% sig line
yMax = max(data(:));
y = yMax + 1;
h = 0.6;
addSigLine(1, 3, y, h, '*');   % Rest vs Motion
addSigLine(2, 3, y + 2.4, h, '***'); % ball vs Motion

% 3. MAR - Short Channels
subplot(2, 2, 3);
data = [final_stats.Rest.MR_Short, final_stats.Task.MR_Short, final_stats.Motion.MR_Short];
boxplot(data, 'Labels', labels, 'Symbol', 'o');
title('MAR - Shorter Channels', 'FontSize', 16, 'FontWeight', 'bold'); 
ylabel('Artifact Ratio (%)', 'FontSize', 16, 'FontWeight', 'bold'); 
ylim([-0.1 3.1]);
grid on;
% sig line
yMax = max(data(:));
y = yMax + 0.15;
h = 0.05;
addSigLine(1, 3, y, h, '***');   % Rest vs Motion
addSigLine(2, 3, y + 0.2, h, '***'); % ball vs Motion

% 4. MAR - Long Channels 
subplot(2, 2, 4);
data = [final_stats.Rest.MR_Long, final_stats.Task.MR_Long, final_stats.Motion.MR_Long];
boxplot(data, 'Labels', labels, 'Symbol', 'o');
title('MAR - Longer Channels', 'FontSize', 16, 'FontWeight', 'bold'); 
ylabel('Artifact Ratio (%)', 'FontSize', 16, 'FontWeight', 'bold'); 
ylim([-0.1 3.1]);
grid on;
% sig line
yMax = max(data(:));
y = yMax + 0.15;
h = 0.05;
addSigLine(1, 3, y, h, '***');   % Rest vs Motion
addSigLine(2, 3, y + 0.2, h, '***'); % ball vs Motion



%% Channel Distance Analysis (Paired T-test)
% 1. Separate all SNR into shorter and longer channel
all_SNR_Short = [final_stats.Rest.SNR_Short];
all_SNR_Long  = [final_stats.Rest.SNR_Long];

% 2. Calculate Differences (Short - Long)
diff_SNR = all_SNR_Short - all_SNR_Long;

% 3. Check Normality of Differences (Crucial Step!)
[h_norm, p_norm] = lillietest(diff_SNR);
fprintf('\n==============================================\n');
fprintf('Paired T-test Normality:\n');
fprintf('p-value = %.4f\n', p_norm)

% 4. Run Paired T-test (Directional: Short > Long)
% 'Tail', 'right' means we test if all_SNR_Short > all_SNR_Long
[~, p_snr, ~, stats_snr] = ttest(all_SNR_Short, all_SNR_Long, 'Tail', 'right');

% Hypothesis Test Shorter SNR > Longer SNR
fprintf('   Mean Shorter: %.2f  |  Mean Longer: %.2f\n', mean(all_SNR_Short), mean(all_SNR_Long));
fprintf('   Mean Diff:    %.2f\n', mean(diff_SNR));
fprintf('   t(%d) = %.2f, p = %.5e (One-tailed)\n', stats_snr.df, stats_snr.tstat, p_snr);
%% Plot Channel Distance Comparison
figure('Color', 'w', 'Position', [150, 150, 800, 500], 'Name', 'Channel Distance Effect');

data_box_snr = [all_SNR_Short, all_SNR_Long]; 

boxplot(data_box_snr, 'Labels', {'Shorter', 'Longer'}, 'Symbol', 'o');
title('SNR Comparison', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('SNR (dB)', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
ylim([20 45]); 

y_lims = ylim;
y_max = max(data_box_snr(:));
y_line = y_max + (y_lims(2)-y_lims(1))*0.05; 
sig_label = '***';
addSigLine(1, 2, y_line, (y_lims(2)-y_lims(1))*0.02, sig_label);

