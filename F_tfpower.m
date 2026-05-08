% Calculate time frequency reconstruction of power from EEG
% amb895 05/07/2026
function [power, times, freqs, varargout] = F_tfpower(EEG,varargin)
% inputs
% EEG = EEG structure from EEGLAB
% optional inputs
% 'channels' = vector of channel numbers
% 'evoked' = 0, not evoked; 1, evoked
% 'itc' = 0, do not compute itc; 1, compute itc
% outputs
% power - freq x time X trials X channels matrix
% times
% freqs
% optional outputs
% itc - freq x time x channels

% defining optional input and output variables
args = reshape(varargin,2,[]);
p = struct(args{:});
nout = max(nargout,1) - 3; % 3 variables that are alway output

% EEG default settings
nchans = EEG.nbchan;
EEGdata = EEG.data;
channels = 1:64;
times = EEG.times;
pnts = EEG.pnts;
ntrls = EEG.trials;
Fs = EEG.srate;
run_itc = 0;
evoked=0;

% define parameters based on inputs
for i = 1:size(args,2)
    param = args{1,i};
    value = args{2,i};
    param = lower(param);
    switch param
        case 'channels' % only run specific channels if specified
            nchans = length(value);
            channels = value;
            EEGdata = EEG.data(value,:,:);
        case 'evoked'
            if value
                evoked =1;
                EEGdata = mean(EEGdata,3); % channels x time
                ntrls = 1;
            else 
                evoked=0;
                EEGdata = EEG.data; % channels x time x trials
                ntrls = EEG.ntrials;
            end
        case 'itc'
            if value
                run_itc = 1;
            else
                run_itc = 0;
            end
    end
end

% Wavelet and convolution parameters
nfreqs = 150;
freq_range = [3 70];
freqs = linspace(freq_range(1), freq_range(2), nfreqs);
cyc = linspace(3,35,nfreqs);
wavet = -2 : 1/Fs : 2; % in seconds
wavet = wavet';
halfw = floor(length(wavet) / 2) + 1;
nconv = pnts + length(wavet) - 1;

%% Manual computation of morlet wavelet reconstruction
if evoked
    EEGdata = permute(EEGdata,[2 1]); % time x channels
else 
    EEGdata = permute(EEGdata, [2 3 1]); % time x trials x channels
end
EEGdataX = fft(EEGdata,nconv,1);

% preallocate temp_power and F
F = zeros(nfreqs,pnts,ntrls,nchans);
power = zeros(nfreqs,pnts,ntrls,nchans);

% loop over frequencies
for fi = 1:nfreqs
    % define wavelet in time domain
    gaus_width = cyc(fi) ./ (2 * pi *freqs(fi));
    wave = exp(2 * 1i * pi * freqs(fi) *wavet) .* exp(-wavet.^2 / (2*gaus_width^2));
    A = 1 / sqrt(2 * pi* freqs(fi)); % normalization coeff
    wave = A*wave;
    % wavelet in frequency domain
    waveX = fft(wave,nconv,1); 
    % Convolution in freq domain
    temp_full = ifft(EEGdataX.*waveX,[],1);
    % trim
    temp = temp_full(halfw:halfw + pnts -1,:,:);
    % save complex reconstruction if running ITC
    if run_itc
        F(fi,:,:,:) = temp;
    end
    % compute power
    temp_power = abs(temp).^2;
    power(fi,:,:,:) = temp_power;
end

% comptue ITC if specified
if run_itc
    itc = squeeze(abs(mean(F ./ abs(F),3)));
    varargout = {itc};
end

end