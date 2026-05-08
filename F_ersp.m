% Calculate ERSP 
% amb895 05/07/2026
function [ERSP_dB, eventtimes] = F_ersp(power,times,baseline_times,event_times,varargin)
% inputs
% power : freqs x times x trials x channels
% times : full epoch times
% baseline_times : [startbase, endbase]
% event_times : [startevent, endevent]
% optional inputs
% base_norm = 'singletrial' baseline normalization or 
%               'avgtrial' average trial normalization for each trial then average across trials
%               'avgtrial_ers' average trial normlization of ERS (default)

% getting optional inputs
args = reshape(varargin,2,[]);
p = struct(args{:});

% default baseline normalization is average trials
base_norm = 3;

for i = 1:size(args,2)
    param = args{1,i};
    value = args{2,i};
    param = lower(param);
    switch param
        case 'base_norm'
            if strcmpi(value,'singletrial')
                base_norm = 1;
            elseif strcmpi(value,'avgtrial')
                base_norm = 2;
            elseif strcmpi(value,'avgtrial_ers')
                base_norm = 3;
            end
    end
end

% find baseline start and end indicies
startbase = find(times >= baseline_times(1), 1, 'first');
endbase = find(times >= baseline_times(2), 1, 'first');

% find event start and end indicies
startevent = find(times >= event_times(1),1,'first');
endevent = find(times >= event_times(2),1,'first');
eventtimes = times(startevent:endevent);

if base_norm == 1 % single trial normalization
    % define baseline
    baseline = mean(power(:,startbase:endbase,:,:),2); % freqs x 1 x trials x channels
    % normalize power with baseline
    ERSP_percent = squeeze(mean(power ./ baseline ,3));  % freqs x times x channels
    % convert to dB
    ERSP_dB = 10*log10(ERSP_percent);
    % just event ersp
    ERSP_dB = ERSP_dB(:,startevent:endevent,:);
elseif base_norm == 2 % average trial normalization w/o ERS
    % define baseline
    baseline = mean(power(:,startbase:endbase,:,:),[2 3]); % freqs x 1 x 1 x channels
    %normalize each trial with average baseline then average across trials
    ERSP_percent = squeeze(mean(power ./ baseline ,3));  % freqs x times x channels
    % convert to dB
    ERSP_dB = 10*log10(ERSP_percent);
    % just event ersp
    ERSP_dB = ERSP_dB(:,startevent:endevent,:);
elseif base_norm == 3
    % define baseline
    baseline = mean(power(:,startbase:endbase,:,:),[2 3]); % freqs x 1 x 1 x channels
    % compute ers
    ERS = mean(power,3); % freqs x times x 1 x channels
    % normalize ERS with baseline
    ERSP_percent = squeeze(ERS ./ baseline);
    % convert to dB
    ERSP_dB = 10*log10(ERSP_percent);
    % just event ersp
    ERSP_dB = ERSP_dB(:,startevent:endevent,:);
end

end