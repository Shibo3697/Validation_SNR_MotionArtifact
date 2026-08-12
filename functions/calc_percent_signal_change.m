function psc = calc_percent_signal_change(d, t, stimOnsets, preWin, postWin)
% CALC_PERCENT_SIGNAL_CHANGE  Percent signal change in light intensity

nStim = length(stimOnsets);
nCh   = size(d, 2);
psc   = zeros(nStim, nCh);

for s = 1:nStim
    onset = stimOnsets(s);

    idxPre  = (t >= onset + preWin(1))  & (t < onset + preWin(2));
    idxPost = (t >= onset + postWin(1)) & (t < onset + postWin(2));

    meanPre  = mean(d(idxPre, :), 1);   % [1 x nCh]
    meanPost = mean(d(idxPost, :), 1);  % [1 x nCh]

    % Percent signal change (Eq. 1 in Yucel et al.)
    psc(s, :) = abs((meanPost - meanPre) ./ meanPre) * 100;
end

end
