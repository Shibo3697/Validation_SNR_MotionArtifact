clear; clc; close all;

%% Configuration
dataFolder  = '/Users/zhoushibo/Desktop/motion task';
sdFilePath  = fullfile(dataFolder, 'sourcedata', 'probe_bilateral.SD');
snirfFiles  = dir(fullfile(dataFolder, '**', '*.snirf'));

distThreshold = 22;       % shorter/longer channel threshold (mm)

taskLevels = ["Rest", "Ball", "Motion"];

%% Channel distance
sd = load(sdFilePath, '-mat');
firstSnirf = SnirfClass(fullfile(snirfFiles(1).folder, snirfFiles(1).name));
ml = firstSnirf.data(1).measurementList;

dists = zeros(1, numel(ml));
for ch = 1:numel(ml)
    srcIdx = ml(ch).sourceIndex;
    detIdx = ml(ch).detectorIndex;
    dists(ch) = norm(sd.SD.SrcPos(srcIdx, :) - sd.SD.DetPos(detIdx, :));
end

idxShorter = dists <= distThreshold;
idxLonger  = dists >  distThreshold;

%% Per-file channel pruning and CV
nFiles = numel(snirfFiles);
subjectID = strings(nFiles, 1);
taskName = strings(nFiles, 1);

retainedShorter = nan(nFiles, 1);
retainedLonger  = nan(nFiles, 1);
cvBeforeShorter = nan(nFiles, 1);
cvBeforeLonger  = nan(nFiles, 1);
cvAfterShorter  = nan(nFiles, 1);
cvAfterLonger   = nan(nFiles, 1);

for k = 1:nFiles
    fprintf('Processing [%d/%d]: %s\n', k, nFiles, snirfFiles(k).name);

    filePath = fullfile(snirfFiles(k).folder, snirfFiles(k).name);
    snirfObj = SnirfClass(filePath);
    rawIntensity = snirfObj.data.dataTimeSeries;

    channelCV = std(rawIntensity, 0, 1, 'omitnan') ./ ...
        mean(rawIntensity, 1, 'omitnan');

    [~, mlActAuto] = evalc("hmrR_PruneChannels(snirfObj.data, snirfObj.probe, [], [], [0 1e7], 10, [0.0 45.0])");
    mlActAuto = mlActAuto{1}(:, 3)';

    retainedShorter(k) = sum(mlActAuto & idxShorter);
    retainedLonger(k)  = sum(mlActAuto & idxLonger);
    cvBeforeShorter(k) = mean(channelCV(idxShorter), 'omitnan');
    cvBeforeLonger(k)  = mean(channelCV(idxLonger),  'omitnan');
    cvAfterShorter(k)  = mean(channelCV(mlActAuto & idxShorter), 'omitnan');
    cvAfterLonger(k)   = mean(channelCV(mlActAuto & idxLonger),  'omitnan');

    subjectID(k) = string(regexp(snirfFiles(k).name, 'sub-[^_]+', 'match', 'once'));
    taskName(k) = get_task_name(snirfFiles(k).name);
end

%% Subject-level
subjects = unique(subjectID, 'stable');
nRows = numel(subjects) * numel(taskLevels);

summaryTask = strings(nRows, 1);
summaryRetainedShorter = nan(nRows, 1);
summaryRetainedLonger  = nan(nRows, 1);
summaryCvBeforeShorter = nan(nRows, 1);
summaryCvBeforeLonger  = nan(nRows, 1);
summaryCvAfterShorter  = nan(nRows, 1);
summaryCvAfterLonger   = nan(nRows, 1);

row = 0;
for s = 1:numel(subjects)
    for t = 1:numel(taskLevels)
        row = row + 1;
        mask = (subjectID == subjects(s)) & (taskName == taskLevels(t));

        summaryTask(row) = taskLevels(t);
        summaryRetainedShorter(row) = mean(retainedShorter(mask), 'omitnan');
        summaryRetainedLonger(row)  = mean(retainedLonger(mask),  'omitnan');
        summaryCvBeforeShorter(row) = mean(cvBeforeShorter(mask), 'omitnan');
        summaryCvBeforeLonger(row)  = mean(cvBeforeLonger(mask),  'omitnan');
        summaryCvAfterShorter(row)  = mean(cvAfterShorter(mask),  'omitnan');
        summaryCvAfterLonger(row)   = mean(cvAfterLonger(mask),   'omitnan');
    end
end

%% Print summary
fprintf('\n--- Channel Pruning Statistics ---\n');
for t = 1:numel(taskLevels)
    mask = taskName == taskLevels(t);
    fprintf('\n%s\n', taskLevels(t));
    print_pruning_total('Shorter', retainedShorter(mask), sum(idxShorter));
    print_pruning_total('Longer ', retainedLonger(mask),  sum(idxLonger));
end

fprintf('\n--- Ball-Squeezing Run-Level c.v. Before Channels Prune ---\n');
fileNames = string({snirfFiles.name})';
for runIdx = 1:3
    runMask = (taskName == "Ball") & contains(fileNames, sprintf('run-%02d', runIdx));
    fprintf('Run %d: Shorter %.4f +/- %.4f | Longer %.4f +/- %.4f\n', ...
        runIdx, ...
        mean(cvBeforeShorter(runMask), 'omitnan'), std(cvBeforeShorter(runMask), 'omitnan'), ...
        mean(cvBeforeLonger(runMask),  'omitnan'), std(cvBeforeLonger(runMask),  'omitnan'));
end

%% Plot
figure('Color','w','Position',[100 100 1200 360]);
tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');

nexttile;
plot_task_pair_box(summaryTask, summaryRetainedShorter, summaryRetainedLonger, ...
    taskLevels);
ylabel('No. of channels with c.v. < 0.1','FontSize',12,'FontWeight','bold');

nexttile;
plot_task_pair_box(summaryTask, summaryCvBeforeShorter, summaryCvBeforeLonger, ...
    taskLevels);
ylabel('c.v. before channels prune','FontSize',12,'FontWeight','bold');

nexttile;
plot_task_pair_box(summaryTask, summaryCvAfterShorter, summaryCvAfterLonger, ...
    taskLevels);
ylabel('c.v. after channels prune','FontSize',12,'FontWeight','bold');
legend({'Shorter','Longer'}, 'Location','northeast');

%% Plot ball-squeezing run-level c.v.
figure('Color','w','Position',[160 160 520 420]);
hold on;
for runIdx = 1:3
    runMask = (taskName == "Ball") & contains(fileNames, sprintf('run-%02d', runIdx));
    boxchart((runIdx - 0.16) * ones(sum(runMask), 1), cvBeforeShorter(runMask), ...
        'MarkerStyle','o', 'MarkerColor',[0.2 0.6 0.9], ...
        'BoxFaceColor', [0.2 0.6 0.9], 'BoxWidth', 0.16);
    boxchart((runIdx + 0.16) * ones(sum(runMask), 1), cvBeforeLonger(runMask), ...
        'MarkerStyle','o', 'MarkerColor',[0.9 0.4 0.3], ...
        'BoxFaceColor', [0.9 0.4 0.3], 'BoxWidth', 0.16);
end
set(gca, 'XTick', 1:3, 'XTickLabel', {'Run 1','Run 2','Run 3'});
ylabel('c.v. before channels prune','FontSize',12,'FontWeight','bold');
legend({'Shorter','Longer'}, 'Location','northeast');
grid on;

%% Local functions
function taskName = get_task_name(fileName)
    if contains(fileName, 'task-Rest', 'IgnoreCase', true)
        taskName = "Rest";
    elseif contains(fileName, 'task-Motion', 'IgnoreCase', true)
        taskName = "Motion";
    else
        taskName = "Ball";
    end
end

function print_pruning_total(labelText, retainedVals, nChannelsPerFile)
    nFiles = sum(~isnan(retainedVals));
    totalChannels = nFiles * nChannelsPerFile;
    retainedChannels = sum(retainedVals, 'omitnan');
    prunedChannels = totalChannels - retainedChannels;
    prunedPercent = prunedChannels / totalChannels * 100;

    fprintf('%s: total=%d | pruned=%d | pruned %.2f%%\n', ...
        labelText, totalChannels, prunedChannels, prunedPercent);
end

function plot_task_pair_box(taskLabels, shorterVals, longerVals, taskLevels)
    hold on;
    for t = 1:numel(taskLevels)
        mask = taskLabels == taskLevels(t);
        boxchart((t - 0.16) * ones(sum(mask), 1), shorterVals(mask), ...
            'MarkerStyle','o', 'MarkerColor',[0.2 0.6 0.9], ...
            'BoxFaceColor', [0.2 0.6 0.9], 'BoxWidth', 0.16);
        boxchart((t + 0.16) * ones(sum(mask), 1), longerVals(mask), ...
            'MarkerStyle','o', 'MarkerColor',[0.9 0.4 0.3], ...
            'BoxFaceColor', [0.9 0.4 0.3], 'BoxWidth', 0.16);
    end
    format_task_axis(taskLevels);
end

function format_task_axis(taskLevels)
    set(gca, 'XTick', 1:numel(taskLevels), 'XTickLabel', taskLevels);
    grid on;
end
