function mask = calc_thres_mask_tfclusters_2d(tfce_scores,perm_tfcescores,p_thres)
% determine pth percentile
percentile = 1-p_thres;

% determine pth percentile for permuted TFCE scores
n = length(perm_tfcescores)*percentile;

% sort permuted TFCE scores in acending order
sorted_permTFCEscores = sort(perm_tfcescores);

% Threshold for 95th percentile
thres = sorted_permTFCEscores(n);

% Threshold Mask
mask = tfce_scores > thres;
end