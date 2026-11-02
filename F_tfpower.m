% Calculate time frequency reconstruction of power from EEG
% amb895 05/07/2026
function [power, times, freqs, varargout] = F_tfpower(EEG, channels, run_evoked, run_itc, run_amplitude)
% inputs
% EEG = EEG structure from EEGLAB
% inputs
%   EEG structure
% optional inputs
%   channels = vector of channel numbers
%   run_evoked = 0, not evoked; 1, evoked
%   run_itc = 0, do not compute itc; 1, compute itc
%   run_amplitude = 0, do not return amplitude envelope; 1, return amplitdue envelope
% outputs
%   power - freq x time X trials X channels 
%   times
%   freqs
% optional outputs
%   itc - freq x time x channels
%   amplitude - freq x time x trials x channels

arguments
    EEG struct
    channels = 1:64
    run_evoked = 0 
    run_itc = 0
    run_amplitude = 0
end
varargout = {};

% EEG default settings
nchans = length(channels);
EEGdata = EEG.data;
times = EEG.times;
pnts = EEG.pnts;
ntrls = EEG.trials;
Fs = EEG.srate;

% for evoked power spectrum, average across trials in time domain
if run_evoked
    EEGdata = mean(EEGdata,3);
    ntrls = 1;
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
% change shape of matrix for computations
if run_evoked
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
    % compute FWHM in time domain
    midp = dsearchn(wavet,0);
    gauswin = exp(-wavet.^2 / (2*gaus_width^2));
    fwhmT(fi) = wavet(midp - 1 + dsearchn(gauswin(midp:end),0.5)) - wavet(dsearchn(gauswin(1:midp),0.5));
    % wavelet in frequency domain
    waveX = fft(wave,nconv,1); 
    % Convolution in freq domain
    temp_full = ifft(EEGdataX.*waveX,[],1);
    % trim
    temp = temp_full(halfw:halfw + pnts -1,:,:);
    % save complex reconstruction if running ITC
    if run_itc || run_amplitude
        F(fi,:,:,:) = temp;
    end
    % compute power
    temp_power = abs(temp).^2;
    power(fi,:,:,:) = temp_power;
end

% comptue ITC if specified
if run_itc
    itc = squeeze(abs(mean(F ./ abs(F),3)));
    varargout{end+1} = itc;
end

% return amplitude envelopes if specified
if run_amplitude
    amp_env = abs(F);
    varargout{end+1} = amp_env;
end


end