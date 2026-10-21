function idx = first_of_seq(ts, v)
% first_of_seqs - find the index where a repeat of a value occurs
% that is, exclude indexes that are a repeat of immedietly previous value
% useful for status channels onsets
%  returns more than `find(ts==v, 'first')`, but fewer than `find(ts==v)`

    hits = find(ts==v);
    % exclude if previous value is also a hit
    % will fail if first hit is 1
    x = arrayfun(@(i) (ts(i-1) ~= v)*i, hits);
    % get the values that are not repeats
    idx = x(x~=0);
end

