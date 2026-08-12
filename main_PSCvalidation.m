%%  Motion Artifact Signal Quality Validation (Yucel et al., 2014)
clear; clc; close all;

%% ======================== Configuration =================================
dataFolder  = '/Users/zhoushibo/Desktop/motion task';
sdFilePath  = fullfile(dataFolder, 'sourcedata', 'probe_bilateral.SD');
snirfFiles  = dir(fullfile(dataFolder, '**', '*.snirf'));
motionFiles = snirfFiles(contains({snirfFiles.name}, 'task-Motion', 'IgnoreCase', true));

preWin  = [-2 0];   % baseline window (s)
postWin = [0 2];    % post-movement window (s)
ppf     = [6 6];    % partial pathlength factors
distThreshold = 22; % shorter/longer channel threshold (mm)
motionLabels = {'Tilting left', 'Tilting right', 'Saying hello', ...
    'Lifting eyebrows', 'Tilting down', 'Tilting up'};
accCzIdx = [1 2 3];       % X/Y/Z accelerometer channels at Cz
axisLabels = {'X', 'Y', 'Z'};
plotFontSize = 12;

set(groot, 'defaultAxesFontSize', plotFontSize, ...
    'defaultTextFontSize', plotFontSize);

%% ======================== Channel Distance ==============================
% Load source/detector positions from the SD file
sdfile = load(sdFilePath, '-mat');
srcPos = sdfile.SD.SrcPos;
detPos = sdfile.SD.DetPos;

% Read channel metadata
tmpSnirf = SnirfClass(fullfile(motionFiles(1).folder, motionFiles(1).name));
ml = tmpSnirf.data(1).measurementList;
numChannels = numel(ml);
wavelengths = tmpSnirf.probe.wavelengths(:)';


SDpairs = zeros(numChannels, 3);
dists   = zeros(1, numChannels);
for i = 1:numChannels
    SDpairs(i, :) = [ml(i).sourceIndex, ml(i).detectorIndex, ml(i).wavelengthIndex];
    dists(i) = norm(srcPos(ml(i).sourceIndex, :) - detPos(ml(i).detectorIndex, :));
end
% SDpair: column 1 = source, 2 = detector, 3 = wavelength

% detect shorter and longer channel pairs
pairDists = dists(SDpairs(:,3)' == 1);
idx_shorter_pair = pairDists <= distThreshold;
idx_longer_pair  = pairDists >  distThreshold;

% 2 wavelengths and shorter/longer channel
idx_shorter_wl1 = (dists <= distThreshold) & (SDpairs(:,3)' == 1);
idx_shorter_wl2 = (dists <= distThreshold) & (SDpairs(:,3)' == 2);
idx_longer_wl1  = (dists >  distThreshold) & (SDpairs(:,3)' == 1);
idx_longer_wl2  = (dists >  distThreshold) & (SDpairs(:,3)' == 2);

% Data Containers
PSC_all = [];    
trialSubj = [];
trialType = [];
HB_stdHbO = [];
HB_stdHbR = [];
HB_dHbO_all = [];
HB_dHbR_all = [];
accChangeCz = [];

%% ======================== Main Loop =====================================
fprintf('Found %d motion artifact induction task files.\n', numel(motionFiles));
for k = 1:numel(motionFiles)
    fprintf('Processing [%d/%d]: %s\n', k, numel(motionFiles), motionFiles(k).name);

    snirfObj = SnirfClass(fullfile(motionFiles(k).folder, motionFiles(k).name));
    d = snirfObj.data.dataTimeSeries;
    t = snirfObj.data.time;

    % extract event
    stimOnsets = [];
    stimTypes  = [];
    for si = 1:numel(snirfObj.stim)
        markerID = str2double(snirfObj.stim(si).name);
        onsets = snirfObj.stim(si).data(:, 1);
        stimOnsets = [stimOnsets; onsets];
        stimTypes  = [stimTypes; repmat(markerID, numel(onsets), 1)];
    end

    % Metric 1: percent signal change
    psc = calc_percent_signal_change(d, t, stimOnsets, preWin, postWin);
    PSC_all  = [PSC_all; psc];
    trialSubj = [trialSubj; repmat(k, size(psc, 1), 1)];
    trialType = [trialType; stimTypes];

    % Metric 2: whole task HbO/HbR std and trial-level HbO/HbR changes
    [sHbO, sHbR, dHbO_m, dHbR_m] = calc_hb_std(snirfObj, stimOnsets, ppf, postWin);
    HB_stdHbO = [HB_stdHbO; sHbO];
    HB_stdHbR = [HB_stdHbR; sHbR];

    HB_dHbO_all = [HB_dHbO_all; dHbO_m]; 
    HB_dHbR_all = [HB_dHbR_all; dHbR_m];

    % Metric 3: Cz accelerometer X/Y/Z changes
    accChangeCz_m = calc_acc_change(snirfObj, stimOnsets, accCzIdx, preWin, postWin);
    accChangeCz = [accChangeCz; accChangeCz_m];
end

%% ======================== Analysis ======================================
subjects = unique(trialSubj);
uniqueTypes = unique(trialType);
typeLabels = motionLabels(uniqueTypes);

%% Metric 1: PSC by Movement Type
psc_s_wl1 = []; psc_l_wl1 = [];
psc_s_wl2 = []; psc_l_wl2 = [];
pscGroup = [];
pscSubj = [];
% _s_=shorter, _l_=longer, _g_=group

for mt = 1:numel(uniqueTypes)
    for s = 1:numel(subjects)
        mask = (trialType == uniqueTypes(mt)) & (trialSubj == subjects(s));
        psc_s_wl1(end+1,1) = mean(mean(PSC_all(mask, idx_shorter_wl1), 2, 'omitnan'), 'omitnan');
        psc_l_wl1(end+1,1) = mean(mean(PSC_all(mask, idx_longer_wl1),  2, 'omitnan'), 'omitnan');
        psc_s_wl2(end+1,1) = mean(mean(PSC_all(mask, idx_shorter_wl2), 2, 'omitnan'), 'omitnan');
        psc_l_wl2(end+1,1) = mean(mean(PSC_all(mask, idx_longer_wl2),  2, 'omitnan'), 'omitnan');
        pscGroup(end+1,1) = mt;
        pscSubj(end+1,1) = s;
    end
end

fprintf('\n--- Subject-Level PSC by Movement Type ---\n');
for mt = 1:numel(uniqueTypes)
    fprintf('%s: Shorter %.2f%%/%d nm, %.2f%%/%d nm | Longer %.2f%%/%d nm, %.2f%%/%d nm\n', ...
        typeLabels{mt}, ...
        median(psc_s_wl1(pscGroup == mt),'omitnan'), wavelengths(1), ...
        median(psc_s_wl2(pscGroup == mt),'omitnan'), wavelengths(2), ...
        median(psc_l_wl1(pscGroup == mt),'omitnan'), wavelengths(1), ...
        median(psc_l_wl2(pscGroup == mt),'omitnan'), wavelengths(2));
end

% % test outlier
% fprintf('\n--- Maximum Subject-Level PSC Data Points ---\n');
% pscVars = {psc_s_wl1, psc_l_wl1, psc_s_wl2, psc_l_wl2};
% pscNames = {'760 nm Shorter', '760 nm Longer', '850 nm Shorter', '850 nm Longer'};
% for i = 1:numel(pscVars)
%     [maxVal, idx] = max(pscVars{i});
%     fileName = motionFiles(subjects(pscSubj(idx))).name;
%     subjectID = regexp(fileName, 'sub-\d+', 'match', 'once');
%     fprintf('%s: %.2f%%, subject %s, motion %s, file %s\n', ...
%         pscNames{i}, maxVal, subjectID, typeLabels{pscGroup(idx)}, fileName);
% end

%% Metric 2a: Whole-Dataset Std of HbO/HbR Changes
sHbO_shorter = mean(HB_stdHbO(:, idx_shorter_pair), 2);
sHbO_longer  = mean(HB_stdHbO(:, idx_longer_pair),  2);
sHbR_shorter = mean(HB_stdHbR(:, idx_shorter_pair), 2);
sHbR_longer  = mean(HB_stdHbR(:, idx_longer_pair),  2);

fprintf('\n--- Subject-Level Whole-Dataset Std of HbO/HbR Changes ---\n');
fprintf('HbO: Shorter %.4f +/- %.4f | Longer %.4f +/- %.4f uM\n', ...
    mean(sHbO_shorter), std(sHbO_shorter), mean(sHbO_longer), std(sHbO_longer));
fprintf('HbR: Shorter %.4f +/- %.4f | Longer %.4f +/- %.4f uM\n', ...
    mean(sHbR_shorter), std(sHbR_shorter), mean(sHbR_longer), std(sHbR_longer));

%% Metric 2b: HbO/HbR Changes by Movement Type
hboS = []; hboL = []; hbrS = []; hbrL = []; hbGroup = [];

for mt = 1:numel(uniqueTypes)
    for s = 1:numel(subjects)
        hmask = (trialType == uniqueTypes(mt)) & (trialSubj == subjects(s));
        dHbO = HB_dHbO_all(hmask, :);
        dHbR = HB_dHbR_all(hmask, :);
        hboS(end+1,1) = mean(mean(dHbO(:, idx_shorter_pair), 2, 'omitnan'), 'omitnan');
        hboL(end+1,1) = mean(mean(dHbO(:, idx_longer_pair),  2, 'omitnan'), 'omitnan');
        hbrS(end+1,1) = mean(mean(dHbR(:, idx_shorter_pair), 2, 'omitnan'), 'omitnan');
        hbrL(end+1,1) = mean(mean(dHbR(:, idx_longer_pair),  2, 'omitnan'), 'omitnan');
        hbGroup(end+1,1) = mt;
    end
end

fprintf('\n--- Subject-Level HbO/HbR Changes by Movement Type ---\n');
for mt = 1:numel(uniqueTypes)
    fprintf('%s: HbO Shorter %.4f uM | HbO Longer %.4f uM | HbR Shorter %.4f uM | HbR Longer %.4f uM\n', ...
        typeLabels{mt}, ...
        median(hboS(hbGroup == mt), 'omitnan'), ...
        median(hboL(hbGroup == mt), 'omitnan'), ...
        median(hbrS(hbGroup == mt), 'omitnan'), ...
        median(hbrL(hbGroup == mt), 'omitnan'));
end

%% Metric 3: Accelerometer (CZ) Changes by Movement Type
czByType = [];
accGroup = [];

for mt = 1:numel(uniqueTypes)
    for s = 1:numel(subjects)
        mask = (trialType == uniqueTypes(mt)) & (trialSubj == subjects(s));
        czByType(end+1, :) = mean(accChangeCz(mask, :), 1, 'omitnan');
        accGroup(end+1, 1) = mt;
    end
end

fprintf('\n--- Subject-Level Cz Accelerometer Change by Movement Type ---\n');
for mt = 1:numel(uniqueTypes)
    vals = median(czByType(accGroup == mt, :), 1, 'omitnan');
    fprintf('%s: %s %.4f | %s %.4f | %s %.4f\n', ...
        typeLabels{mt}, axisLabels{1}, vals(1), axisLabels{2}, vals(2), axisLabels{3}, vals(3));
end

%% Metric 4: Spearman Correlation Between Normalized Accelerometer Change and Log-Transformed Optical PSC
pscByWavelengthDistance = {psc_s_wl1, psc_l_wl1; psc_s_wl2, psc_l_wl2};
distanceLabels = {'Shorter', 'Longer'};

rhoCorr = nan(numel(wavelengths), numel(distanceLabels), numel(axisLabels));% Spearman rho
pCorr = nan(numel(wavelengths), numel(distanceLabels), numel(axisLabels));% pvalue
pCorrBonf = nan(numel(wavelengths), numel(distanceLabels), numel(axisLabels));% Bonferroni-adjusted pvalue
nCorr = nan(numel(wavelengths), numel(distanceLabels), numel(axisLabels));% num of points used
nCorrTests = numel(pCorr);

fprintf('\n--- Subject-Level Spearman Correlation');
for wl = 1:numel(wavelengths)
    fprintf('\n%d nm PSC\n', wavelengths(wl));
    for distIdx = 1:numel(distanceLabels)
        fprintf('%s:', distanceLabels{distIdx});
        for ax = 1:numel(axisLabels)
            pscVals = pscByWavelengthDistance{wl, distIdx};
            validMask = ~isnan(czByType(:, ax)) & ~isnan(pscVals);
            nCorr(wl, distIdx, ax) = sum(validMask);% record how many n are used
            xNorm = normalize(log10(czByType(validMask, ax)), 'zscore');
            yNorm = normalize(log10(pscVals(validMask)), 'zscore');
            [rhoCorr(wl, distIdx, ax), pCorr(wl, distIdx, ax)] = ...
                corr(xNorm, yNorm, 'Type','Spearman');
            pCorrBonf(wl, distIdx, ax) = min(pCorr(wl, distIdx, ax) * nCorrTests, 1);
            fprintf(' Cz-%s rho=%.3f, p=%.4g, p_Bonf=%.4g, n=%d', ...
                axisLabels{ax}, rhoCorr(wl, distIdx, ax), pCorr(wl, distIdx, ax), ...
                pCorrBonf(wl, distIdx, ax), ...
                nCorr(wl, distIdx, ax));
            if ax < numel(axisLabels)
                fprintf(' |');
            end
        end
        fprintf('\n');
    end
end

%% ======================== Plots =========================================
figure('Color','w','Position',[80 80 1400 900]);
tiledlayout(2, 2, 'TileSpacing','compact','Padding','compact');

% 1a: PSC by movement type.
nexttile;
plot_paired_box(pscGroup, psc_s_wl1, psc_l_wl1, typeLabels, ...
    sprintf('PSC at %d nm (%%)', wavelengths(1)), plotFontSize);

nexttile;
plot_paired_box(pscGroup, psc_s_wl2, psc_l_wl2, typeLabels, ...
    sprintf('PSC at %d nm (%%)', wavelengths(2)), plotFontSize);

% 1b: trial-window HbO/HbR changes by movement type.
nexttile;
plot_paired_box(hbGroup, hboS, hboL, typeLabels, 'Change in HbO (\muM)', plotFontSize);
yl = ylim;
ylim([yl(1), max(yl(2), 10)]);

nexttile;
plot_paired_box(hbGroup, hbrS, hbrL, typeLabels, 'Change in HbR (\muM)', plotFontSize);
yl = ylim;
ylim([yl(1), max(yl(2), 10)]);

% 2. whole-dataset standard deviation of HbO/HbR changes
figure('Color','w','Position',[120 120 520 420]);
barMeans = [mean(sHbO_shorter), mean(sHbO_longer); ...
    mean(sHbR_shorter), mean(sHbR_longer)];
barSD = [std(sHbO_shorter), std(sHbO_longer); ...
    std(sHbR_shorter), std(sHbR_longer)];
b = bar(barMeans, 'grouped');
b(1).FaceColor = [0.2 0.6 0.9];
b(2).FaceColor = [0.9 0.4 0.3];
hold on;
for i = 1:numel(b)
    errorbar(b(i).XEndPoints, barMeans(:, i), barSD(:, i), ...
        'k','LineStyle','none','LineWidth',1.2);
end
ylabel('Standard Deviation (\muM)','FontSize',plotFontSize,'FontWeight','bold');
set(gca,'XTickLabel',{'HbO','HbR'});
legend({'Shorter','Longer'},'Location','northwest','FontSize',plotFontSize);
set(gca,'FontSize',plotFontSize);
grid on;

% 3. accelerometer changes by movement type
figure('Color','w','Position',[100 100 1500 430]);
boxWidth = 0.16;
offset = 0.16;
boxchart(accGroup - offset, czByType(:, 1), 'MarkerStyle','o', ...
    'BoxFaceColor','r','BoxWidth',boxWidth,'BoxFaceAlpha',0.25);
hold on;
boxchart(accGroup, czByType(:, 2), 'MarkerStyle','o', ...
    'BoxFaceColor','g','BoxWidth',boxWidth,'BoxFaceAlpha',0.25);
boxchart(accGroup + offset, czByType(:, 3), 'MarkerStyle','o', ...
    'BoxFaceColor','b','BoxWidth',boxWidth,'BoxFaceAlpha',0.25);
ylabel('Cz Acceleration Change (m/s^2)','FontSize',plotFontSize,'FontWeight','bold');
legend({'X','Y','Z'},'Location','northwest','FontSize',plotFontSize);
set(gca,'XTick',1:numel(typeLabels),'XTickLabel',typeLabels);
xlim([0.5 numel(typeLabels)+0.5]);
xtickangle(30);
set(gca,'FontSize',plotFontSize);
grid on;

% 4. correlation between normalized Cz accelerometer changes and log-transformed optical PSC
figure('Color','w','Position',[80 80 1350 950]);
tCorr = tiledlayout(numel(wavelengths) * numel(distanceLabels), numel(axisLabels), ...
    'TileSpacing','compact','Padding','compact');

for wl = 1:numel(wavelengths)
    for distIdx = 1:numel(distanceLabels)
        rowIdx = (wl - 1) * numel(distanceLabels) + distIdx;
        for ax = 1:numel(axisLabels)
            nexttile((rowIdx - 1) * numel(axisLabels) + ax);
            pscVals = pscByWavelengthDistance{wl, distIdx};
            validMask = ~isnan(czByType(:, ax)) & ~isnan(pscVals);
            xNorm = normalize(log10(czByType(validMask, ax)), 'zscore');
            yNorm = normalize(log10(pscVals(validMask)), 'zscore');
            scatter(xNorm, yNorm, 50, ...
                accGroup(validMask), 'filled', 'MarkerFaceAlpha', 0.75);
            if rowIdx == 1
                title(sprintf('Cz-%s: ρ=%.3f, p=%.4g, n=%d', ...
                    axisLabels{ax}, rhoCorr(wl, distIdx, ax), ...
                    pCorrBonf(wl, distIdx, ax), nCorr(wl, distIdx, ax)), ...
                    'FontSize', plotFontSize, 'FontWeight', 'bold');
            else
                title(sprintf('ρ=%.3f, p=%.4g, n=%d', ...
                    rhoCorr(wl, distIdx, ax), pCorrBonf(wl, distIdx, ax), ...
                    nCorr(wl, distIdx, ax)), ...
                    'FontSize', plotFontSize, 'FontWeight', 'bold');
            end
            if ax == 1
                ylabel(sprintf('%s %d nm', distanceLabels{distIdx}, wavelengths(wl)), ...
                    'FontSize', plotFontSize, 'FontWeight', 'bold');
            end
            if rowIdx == numel(wavelengths) * numel(distanceLabels)
                xlabel(sprintf('Cz-%s', axisLabels{ax}), ...
                    'FontSize', plotFontSize, 'FontWeight', 'bold');
            end
            yl = ylim;
            ylim([yl(1), 3]);
            clim([0.5 numel(typeLabels)+0.5]);
            set(gca,'FontSize',plotFontSize);
            grid on;
            box off;
        end
    end
end

xlabel(tCorr, 'Normalized log_{10}(Cz acceleration change)', 'FontSize', plotFontSize, 'FontWeight', 'bold');
ylabel(tCorr, 'Normalized log_{10}(PSC)', 'FontSize', plotFontSize, 'FontWeight', 'bold');
colormap(lines(numel(typeLabels)));
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Ticks = 1:numel(typeLabels);
cb.TickLabels = typeLabels;
cb.Label.String = 'Motion type';
cb.FontSize = plotFontSize;
cb.Label.FontSize = plotFontSize;
cb.Label.FontWeight = 'bold';

%% ======================== Local Functions ===============================
function plot_paired_box(group, shorterVals, longerVals, labels, yLabelText, fontSize)
    boxchart(group - 0.18, shorterVals, 'MarkerStyle','o', ...
        'BoxFaceColor',[0.2 0.6 0.9], 'BoxWidth',0.16);
    hold on;
    boxchart(group + 0.18, longerVals, 'MarkerStyle','o', ...
        'BoxFaceColor',[0.9 0.4 0.3], 'BoxWidth',0.16);
    ylabel(yLabelText,'FontSize',fontSize,'FontWeight','bold');
    legend({'Shorter','Longer'},'Location','northwest','FontSize',fontSize);
    set(gca,'XTick',1:numel(labels),'XTickLabel',labels);
    xlim([0.5 numel(labels)+0.5]);
    xtickangle(30);
    set(gca, 'YScale', 'log');
    set(gca,'FontSize',fontSize);
    grid on;
end
