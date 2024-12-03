function [] = preprocessing_pipelineAS(inputfile, outpath, lowBP, topBP, FLAG, condition, task, varargin)
% Making changes to original preprocessing_pipeline.m for AS data
% Bad channel rejection in original preprocessing_pipeline.m removed events

% what file are we using
if ~exist(inputfile,'file'), error('inputfile "%s" does not exist!', inputfile), end
[d, currentName, ext ] = fileparts(inputfile);
parts = split(currentName,'_');
subid = str2double(parts{1});
scandate = str2double(parts{2});

%% cap locations
eeglabpath = fileparts(which('eeglab'));
cap_location = fullfile(eeglabpath,'/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
if ~exist(cap_location, 'file'), error('cannot find file for 128 channel cap: %s', cap_location), end
correction_cap_location = hera('Projects/7TBrainMech/scripts/eeg/Shane/resources/ChanLocMod128to64.ced');
if ~exist(correction_cap_location, 'file'), error('cannot find file for correction 128 channel cap: %s', correction_cap_location), end

%% Files
subj_files = file_locs(inputfile, outpath, task);

% to know how far your script is with running
fprintf('==========\n%s:\n\t Initial Preprocessing(%s,%f,%f,%s,%s)\n',...
    currentName, inputfile, lowBP, topBP, outpath,task)

% where to save things
filter_folder = 'filtered';
chanrj_folder = 'channels_rejected';
epoch_folder = 'epoched';
icawholein_folder ='rerefwhole';
epoch_rj_marked_folder = 'marked_epochs';
epochrj_folder = 'rejected_epochs';
icaout = fullfile(outpath, 'ICA');
icawholeout = fullfile(outpath, 'ICAwhole');
icawholeouthomog = fullfile(outpath, 'AfterWhole/ICAwholeClean_homogenize');

% and what files will we create
rerefwhole_name = [currentName '_rerefwhole'];
chrm_name   = [currentName '_badchannelrj'];
epochrj_name = [currentName '_epochs_rj'];
epochrj = fullfile(outpath, epochrj_folder, [epochrj_name '.set']); % not sure what this is saving for
% FIXME: these is not actually used!? but is recoreded in data_removed WF20190911
datarm_name = [currentName '_baddatarj'];
epochrm_name = [currentName '_badepochrj'];
icawholeout_name = [currentName '_rerefwhole_ICA'];

% channel names
commonPlus = {'AFz','C1','C2','C3','C4','C5','C6','CP1','CP2','CP3','CP4',...
    'CP5','CP6','CPz','Cz','F1','F2','F3','F4','F5','F6','F7','F8','FC1',...
    'FC2','FC3','FC4','FC5','FC6','FCz','Fp1','Fp2','FT10','FT9','Fz','I1',...
    'I2','O1','O2','Oz','P1','P10','P2','P3','P4','P5','P6','P7','P8','P9',...
    'PO10','PO3','PO4','PO7','PO8','PO9','POz','Pz','T7','T8',...
    'AF8','AF7','AF4','AF3'};

% checking if whole ICA run already exists (subject already preprocessed)
if condition == 1
    icawholeoutFile = fullfile(icawholeout, [icawholeout_name '.set']);
else
    icawholeoutFile = 'no';
end

if exist(icawholeoutFile, 'file')
    warning('%s already complete (have "%s")! todo load from file', currentName, icawholeout_name)
    return
end

% loading filtered EEG if already created
if condition == 1 
    xEEG = load_if_exists(subj_files.filter);
end

% making filtered EEG current EEG if already did filtering
if isstruct(xEEG)
    [ALLEEG EEG CURRENTSET] = pop_newset([], xEEG, 0,...
        'setname',currentName,...
        'gui','off');
else % if filtered EEG does not exist running rereferencing, filtering, and resampling
    %% load EEG set and re-reference to mastoids
     EEG = pop_loadset(inputfile);

    if size(EEG.data,1) < 100
        % [65 66] are mastoid externals
        Flag128 = 0;
        EEG = pop_reref(EEG, [65 66]);
        EEG = eeg_checkset(EEG);
    else
        Flag128 = 1;
        %[129 130] are the mastoid externals for the 128 electrode
        EEG = pop_reref(EEG, [129 130]);
        EEG = eeg_checkset(EEG);
    end
    
    %stores EEG set in ALLEEG, give setname
    ALLEEG = [];
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,...
        'setname',currentName,...
        'gui','on');
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    
    EEG.subject = currentName(1:findstr(currentName,'anti')-2);
    EEG.condition =  currentName(findstr(currentName,'anti'):end);
    
    %% Filtering
    
    EEG = pop_eegfiltnew(EEG, lowBP, topBP, 3380, 0, [], 0);
    % inputs: (EEG, locutoff, hicutoff, filtorder, revfilt, usefft, plotfreqz, minphase);
    % revfilt- reverse filter 0=bandpass, 1=notch
    % usefft- [] ignore (use fft to filter)
    % plotfreqz- plot frequency bode plots for filter 0=don't plot, 1=plot
    % filtorder = 3380 - filter order (filter length - 1). Mandatory even. performing 3381 point bandpass filtering.
    
    %give a new setname and overwrite unfiltered data
    EEG = pop_editset(EEG,'setname',[currentName '_bandpass_filtered']);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    
    %% Rename EEG.events.types to correct, incorrect, or dropped
    % remove incorrect and dropped trials from EEG.data
    addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/eog_cal') % need to load AS score data table
    load('eeg_anti.mat','datatable')
    subsetdatatable = datatable(:,{'LunaID','ScanDate','Trial','Correct','trialOnsetIndex'});
    % subsetting table so it is just the subject and scan date of interest
    subIdx = find(subsetdatatable{:,'LunaID'}==subid);
    scanIdx = find(subsetdatatable{subIdx,'ScanDate'}==scandate);
    subscanIdx = subIdx(scanIdx);
    subscanTrialTable = subsetdatatable(subscanIdx,:);
    trialtable = subscanTrialTable{:,3:5}; % col1=trial#, col2=score, col3=trialonsetindex
    ntrials = length(trialtable);
    
    % making all event types strings
    for nevents=1:length({EEG.event(:).type})
        if ischar(EEG.event(nevents).type)
            continue
        else
            EEG.event(nevents).type = num2str(EEG.event(nevents).type);
        end
    end
    
    % renaming events if they were correct, incorrect, or a dropped trial
    for triali = 1:ntrials
        trialOnsetIdx = trialtable(triali,3);
        if triali ~= 40
            trialOffsetIdx = trialtable(triali+1,3)-1;
        elseif triali == 40
            trialOffsetIdx = length(EEG.data);
        end
        eventIdx = find(cell2mat({EEG.event(:).latency}) >= trialOnsetIdx & cell2mat({EEG.event(:).latency}) <= trialOffsetIdx);
        trialEvents = EEG.event(eventIdx);
        trialScore = trialtable(triali,2);
        prepEventIdxTrial = find(cell2mat({trialEvents.type})=='2');
        prepEventIdxWhole = eventIdx(prepEventIdxTrial);

        if trialScore == 1
            % eeg.event.type 2_cor
            EEG.event(prepEventIdxWhole).type = '2_cor';
        elseif trialScore == 0
            % eeg.event.type 2_incor
            EEG.event(prepEventIdxWhole).type = '2_incor';
        elseif isnan(trialScore)
            % eeg.event.type 2_drop
            EEG.event(prepEventIdxWhole).type = '2_drop';
        end
        trialOnsetTimes(triali) = EEG.times(trialOnsetIdx);
        trialOffsetTimes(triali) = EEG.times(trialOffsetIdx);
    end
    
    % remove dropped and incorrect events and data
    EEG_corTrials = pop_rmdat(EEG,{'2_incor','2_drop'},[-0.5 0.5],1);
    boundary_inds = find(strcmp({EEG_corTrials.event.type}, 'boundary'));
    EEG_corTrials.event(boundary_inds) = []; 
    EEG = EEG_corTrials;
    
    %% Resample Data
    
    % Downsample the data to 150Hz using anti-aliasing filter
    EEG = pop_resample(EEG, 150, 0.8, 0.4); 
    %0.8 is fc and 0.4 is df. Default is .9 and .2. We dont know why Alethia changed them
    % df = anti-aliasing filter transition band width
    % fc = anti-aliasing filter cutoff
    EEG = eeg_checkset(EEG);
    
    %change setname
    EEG = pop_editset(EEG,'setname',[currentName '_filtered']);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    %save filtered data
    EEG = pop_saveset( EEG, 'filename',[currentName '_filtered'], ...
        'filepath',fullfile(outpath,filter_folder));
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
end

%% Channels
% remove external channels
EEG = pop_select( EEG,'nochannel',{'EX3' 'EX4' 'EX5' 'EX6' 'EX7' 'EX8' 'EXG1' 'EXG2' 'EXG3' 'EXG4' 'EXG5' 'EXG6' 'EXG7' 'EXG8' 'GSR1' 'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp' 'FT7' 'FT8' 'TP7' 'TP8' 'TP9' 'TP10'});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%import channel locations
EEG=pop_chanedit(EEG, 'lookup', cap_location);

% fixing if 128 channel setup
if size(EEG.data,1) > 100
    EEG = pop_select( EEG,'channel',commonPlus);
    EEG = pop_chanedit(EEG, 'load', {correction_cap_location 'filetype' 'autodetect'});
    % 128    'AF8' --> 64    'AF6'
    % 128    'AF7' --> 64    'AF5'
    % 128    'AF4' --> 64    'AF2'
    % 128    'AF3' --> 64    'AF1'
end
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Bad Channel Rejection
% bad channels are in general the channels that you can't save even if you reject 5-10% of your datapoints
% from eeglab wiki: plot channel spectrum and ID outliers with spectopo()
% 3 options for channel rejection:
% 1. look at standard deviation in bar plots and remove the channels
% 2. kurtosis
% 3. clean_rawdata - this method removed events from EEG set and lost individual anti trials
% maybe use clean_rawdata but need to tweak from original pipeline

% loading rejected channels if already done
if condition == 1
    xEEG = load_if_exists(subj_files.chanrj);
end
% originalEEG has channel locations and external channels removed- need for determining which channels have been rejected
originalEEG = EEG; 

if isstruct(xEEG) % if channels have already been rejected, setting that set as current EEG
    [ALLEEG EEG] = eeg_store(ALLEEG, xEEG, CURRENTSET);
else
    EEG = clean_rawdata(EEG, 8, -1, 0.7, 5, 15, 0.3);
    % clean_rawdata inputs:
    % (EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst,arg_window)
    % see clean_rawdata for description of what inputs mean
    % arg_highpass = -1 : already highpassed data do not need to do again
    % LOOK AT OTHER PARAMETERS
    
    %change setname
    EEG = pop_editset(EEG,'setname', chrm_name);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG);

    EEG = pop_saveset( EEG,'filename', chrm_name, ...
        'filepath',sprintf('%s/%s/',outpath,chanrj_folder));
end

if condition == 1 % loading average re-referenced EEG set if already created
    xEEG = load_if_exists(subj_files.rerefwhole_name);
end

if isstruct(xEEG) % if all ready re-referenced don't do again
    [ALLEEG EEG] = eeg_store(ALLEEG, xEEG, CURRENTSET);
else
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    % what is clean_channel_mask???
    if ~any(find(cellfun (@any,regexpi (fieldnames(EEG.etc), 'clean_channel_mask'))))
        EEG.etc.clean_channel_mask=42;
    else
    end
    
    %save the channels that were rejected in a variable
    channels_removed{1} = chrm_name; %setname
    channels_removed{2} = setdiff({originalEEG.chanlocs.labels},{EEG.chanlocs.labels}, 'stable');
    channels_removed{3} = find(EEG.etc.clean_channel_mask==0);
    %also save the channels that were rejected in the EEG struct
    EEG.channels_rj = channels_removed{2};
    EEG.channels_rj_nr = length(EEG.channels_rj);
    
% since i disabled clean_window is that why there is no EEG.etc.clean_sample_mask??
    %save the proportion of the dataset that were rejected in a variable
    data_removed{1} = datarm_name; %setname
    data_removed{2} = length(find(EEG.etc.clean_sample_mask==0))/EEG.srate;%
    data_removed{3} = length(find(EEG.etc.clean_sample_mask==0))/length(EEG.etc.clean_sample_mask);%
    %also save the data that were rejected in the EEG struct
    EEG.data_rj    = data_removed{2};
    EEG.data_rj_nr = data_removed{3};
    
    %% Interpolate Channels
    % POSSIBLE PROBLEMS
    %  - injecting extra channels ontop of expected 64 (n>64)
    %  - 128 missing expected labels, adding too few back (n<64)
    if Flag128 == 1 % no clue what is going on here using from original pipeline
        nchan = 64;
        ngood = length(EEG.chanlocs);
        %  128 cap doesn't have exactly the same postions as 64
        % remove 4 that are in the wrong place and reinterpret
        % AND interp any bad channels
        % do this by removing the 4 128weirdos
        % from the already trimmed (no bad channels) in EEG.chanlocs

        need_128interp = [2  3  35  36 ];
        % get the names of those to remove
        n128name = {originalEEG.chanlocs(need_128interp).labels};
        % should always be {'AF5','AF1','AF2','AF6'} ??

        % find where they are in current EEG files (if they haven't already been removed)
        n128here_idx = find(ismember({EEG.chanlocs.labels},n128name));
        % keep those that aren't the ones we matched
        % remove from chanlocs, data and update nbcan
        % WARNING -- who knows what else we should have changed to update the set info!
        keep_idx = setdiff(1:ngood, n128here_idx);
        EEG.chanlocs = EEG.chanlocs(keep_idx);
        EEG.data = EEG.data(keep_idx,:);
        EEG.nbchan = length(keep_idx);

        %EEG_i = pop_interp(EEG, interp_ch, 'spherical');
        fprintf('%d channels in orig; want to interpolate %d bad and move %d\n',...
            originalEEG.nbchan, nchan - ngood, length(need_128interp))
        EEG_i = pop_interp(EEG, originalEEG.chanlocs, 'spherical');

        % could swap these channels (they're close, but not the same)
        % BUT WE DONT
        % 128    'AF7' --> 64    'AF5' In this point channel 2
        % 128    'AF3' --> 64    'AF1' In this point channel 3
        % 128    'AF4' --> 64    'AF2' In this point channel 35
        % 128    'AF8' --> 64    'AF6' In this point channel 36
        % lines above modify channel information and pocition in data to make
        %  it the same for 64 and 128 cap

        % need to do destructive swapping. need a copy
        EEG = EEG_i;
        EEG.chanlocs(2) = EEG_i.chanlocs(3);%EEG.chanlocs(2) must by 'AF1' in 64 cap
        EEG.chanlocs(3) = EEG_i.chanlocs(2);%EEG.chanlocs(3) must by 'AF5' in 64 cap
        %     EEG.chaninfo.filecontent(4,:) = EEG_i.chaninfo.filecontent(3,:); This
        %     is not necessary i think, but just in case...
        %     EEG.chaninfo.filecontent(4,1) = '3';
        %     EEG.chaninfo.filecontent(3,:) = EEG_i.chaninfo.filecontent(4,:);
        %     EEG.chaninfo.filecontent(3,1) = '2';
        EEG.data(2,:) = EEG_i.data(3,:);% ALERT ALERT Lines latelly added
        EEG.data(3,:) = EEG_i.data(2,:);% ATERT ALERT Lines latelly added
    else
        EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
    end
    
    %% Re-reference: Average Reference
    % not sure what FLAG is
    if FLAG
        EEG = pop_reref( EEG, []);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,...
            'setname',[currentName '_avref'],...
            'gui','off');
    end
    
    %save whole rereferenced data for ICA whole
    %save epochs rejected EEG data
    EEG = pop_editset(EEG, 'setname', rerefwhole_name);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    EEG = pop_saveset(EEG, 'filename', rerefwhole_name, 'filepath', fullfile(outpath,icawholein_folder));
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
end
    
%% Whole ICA run
% what is the difference betweeen this and manual ICA?

icawholein = fullfile(outpath, icawholein_folder, [rerefwhole_name '.set']);

if ~exist(subj_files.icawhole, 'file')
    runICAs(icawholein,icawholeout,task)
else
    fprintf('have %s, not rerunning\n', subj_files.icawholeout)
end

end