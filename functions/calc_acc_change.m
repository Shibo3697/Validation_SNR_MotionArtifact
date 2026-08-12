function accChange = calc_acc_change(snirfObj, stimOnsets, accIdx, preWin, postWin)
% CALC_ACC_CHANGE  X/Y/Z accelerometer change for each motion event.
%   Change is absolute difference between post-movement and baseline means.

tAcc = snirfObj.aux(accIdx(1)).time;
accData = [snirfObj.aux(accIdx(1)).dataTimeSeries, ...
    snirfObj.aux(accIdx(2)).dataTimeSeries, ...
    snirfObj.aux(accIdx(3)).dataTimeSeries];

nStim = length(stimOnsets);
accChange = zeros(nStim, 3);

for s = 1:nStim
    onset = stimOnsets(s);

    idxPre = (tAcc >= onset + preWin(1)) & (tAcc < onset + preWin(2));
    idxPost = (tAcc >= onset + postWin(1)) & (tAcc < onset + postWin(2));

    meanPre = mean(accData(idxPre, :), 1, 'omitnan');
    meanPost = mean(accData(idxPost, :), 1, 'omitnan');
    accChange(s, :) = abs(meanPost - meanPre);
end

end
