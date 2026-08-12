function [stdHbO, stdHbR, dHbO_motion, dHbR_motion] = calc_hb_std(snirfObj, stimOnsets, ppf, postWin)
% CALC_HB_STD  HbO/HbR standard deviation and movement-related changes.
%   Uses Homer3 functions:
%       hmrR_Intensity2OD  - intensity to optical density (OD)
%       hmrR_OD2Conc       - OD to HbO/HbR

% --- Step 1: Intensity -> OD (Homer3) ---
dod = hmrR_Intensity2OD(snirfObj.data);

% --- Step 2: OD -> HbO/HbR concentrations  ---
dc = hmrR_OD2Conc(dod, snirfObj.probe, ppf);

% Extract HbO and HbR
dcData = dc.dataTimeSeries;
nCols  = size(dcData, 2);
nPairs = nCols / 3;

HbO = dcData(:, 1:3:end) * 1e6;  % columns 1, 4, 7, ... are HbO, converted to uM
HbR = dcData(:, 2:3:end) * 1e6;  % columns 2, 5, 8, ... are HbR, converted to uM

t = snirfObj.data.time;

% --- Step 3: Std of HbO/HbR over entire dataset ---
stdHbO = std(HbO, 0, 1);  % [1 x nPairs]
stdHbR = std(HbR, 0, 1);  % [1 x nPairs]

% --- Step 4: Mean change in HbO/HbR during each movement ---
nStim = length(stimOnsets);
dHbO_motion = zeros(nStim, nPairs);
dHbR_motion = zeros(nStim, nPairs);

for s = 1:nStim
    onset = stimOnsets(s);
    idxPost = (t >= onset + postWin(1)) & (t < onset + postWin(2));

    % Absolute mean change during motion window
    dHbO_motion(s, :) = abs(mean(HbO(idxPost, :), 1));
    dHbR_motion(s, :) = abs(mean(HbR(idxPost, :), 1));
end

end
