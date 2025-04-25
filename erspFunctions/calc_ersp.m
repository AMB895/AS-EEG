function [suberspdata,subitcdata,subpowbase,subsignalchange,times,freqs,failcode,subID,scanDate,age,sub_ersp_data_name,sub_itc_data_name] = calc_ersp(epoch_path,sub_epoch_filename,sub_ersp_path,subMerge7t,channelTemplate)
% failcodes:
%    1 = not in merge7t
%    2 = zero epochs
%    3 = no age
%    4 = channel labels do not match template

% defining input EEG file
inputfile = [epoch_path,sub_epoch_filename];
[~,currentName,~] = fileparts(inputfile);
splitCurrentName = split(currentName,'_');
subID = cell2mat(splitCurrentName(1));
scanDate = cell2mat(splitCurrentName(2));
lunaid = [subID,'_',scanDate];
sub_ersp_data_name = sprintf('%sersp.mat',lunaid);
sub_itc_data_name = sprintf('%sitc.mat',lunaid);

% determine if subject's ersp has already been calculated
if exist(fullfile(sub_ersp_path,sub_ersp_data_name),'file')
    fprintf('computed ersp for %s\nLoading suberspdata\n',lunaid)
    subStructERSP = load(fullfile(sub_ersp_path,sub_ersp_data_name));
    suberspdata = subStructERSP.suberspdata;
    subpowbase = subStructERSP.subpowbase;
    subsignalchange = subStructERSP.subsignalchange;
    times = subStructERSP.times;
    freqs = subStructERSP.freqs;
    subStructITC = load(fullfile(sub_ersp_path,sub_itc_data_name),'subitcdata');
    subitcdata = subStructITC.subitcdata;
    EEG = pop_loadset(inputfile);
    age = EEG.age;
    failcode = NaN;
    return
else
    % skipping if subject is not in merge 7t
    if isempty(subMerge7t)
        fprintf('\nSubject %s not in Merge7t\n',lunaid)
        failcode = 1;
        suberspdata = [];
        subpowbase = [];
        subitcdata = [];
        subsignalchange=[];
        times = [];
        freqs = [];
        age = [];
        return
    end
    
    % loading subject's epoched EEG
    EEG = pop_loadset(inputfile);
    chanNames = string({EEG.chanlocs.labels}');
    numchans = length(chanNames);
    
    % skipping if subject does not have age
    if isfield(EEG,'age')
        age = EEG.age;
    else
        fprintf('\nSubject %s does not have age\n',lunaid)
        failcode = 3;
        suberspdata = [];
        subpowbase = [];
        subitcdata = [];
        subsignalchange=[];
        times = [];
        freqs = [];
        age=[];
        return
    end
    
    % check that channel labels are the same
    if all(channelTemplate~=chanNames)
        fprintf('\nCheck Channel Labels %s\n',lunaid)
        failcode = 4;
        suberspdata = [];
        subpowbase = [];
        subitcdata = [];
        subsignalchange=[];
        times = [];
        freqs = [];
        return
    end
    
    % skip if subject has zero epochs
    numEpochs = length(EEG.epoch);
    if numEpochs < 1
        fprintf('\nSubject %s has %d epochs\n',lunaid,numEpochs)
        failcode = 2;
        suberspdata = [];
        subpowbase = [];
        subitcdata = [];
        subsignalchange=[];
        times = [];
        freqs = [];
        return
    end
    
    % calculate subject's ersp and itc for each channel
    for currentChan = 1:numchans
        % printing progress
        fprintf('\n\nChannel %d out of %d\n\n',currentChan,numchans)
        % Inputs to newtimef():
        % Required:
        %   data: 3D array (1,frames,trials)
        %   frames: frames/trial; ignored for 3D data
        %   tlimits: time limits of data epochs, NOT subwindow to extract from epcochs (in ms)
        %   Fs: data sampling rate
        %   varwin: [# of wavelet cycles @ highest freq, window size increase]
        % Optional:
        %   'freqs': frequency limits
        %   'baseline': spectral baseline range [-700 -400] from Kai
        %   'plotersp' & 'plotitc' & 'plotphasesign': plotting options
        % Filtered original data with bandpass filter 0.5-70 Hz and notch filter at 60 Hz-> should not look at frequencies above 60 Hz
        [suberspdata(currentChan,:,:),subitcdata(currentChan,:,:),subpowbase(currentChan,:),times,freqs]=newtimef(EEG.data(currentChan,:,:),...
             EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate, [3 0.5], 'freqs', [3 59],...
            'baseline',[-700 -400],'nfreqs',100,'timesout',250,'freqscale','log',...
            'plotersp','off','plotitc','off','plotphasesign','off');
    end
    % subsetting ersp and itc data for just prep period
    time0 = find(times==0);
    time500 = find(times==500);
    suberspdata = suberspdata(:,:,time0:time500);
    subitcdata = subitcdata(:,:,time0:time500);
    times = times(time0:time500);
    failcode = NaN;
    
    % Calculating percent signal change from baseline power
    numTimes = length(times);
    for currentChan = 1:numchans
        chanersp = squeeze(suberspdata(currentChan,:,:));
        chanpowbase = squeeze(subpowbase(currentChan,:))';
        for currentT = 1:numTimes
            chanTersp = squeeze(chanersp(:,currentT));
            subsignalchange(currentChan,:,currentT) = (chanTersp./chanpowbase).*100;
        end
    end
end