function [revisar] = epochAS(inputfile, epochfolder, task, taskdirectory,condition,lowBP,topBP)
revisar = {};

% what file are we using
if ~exist(inputfile,'file'), error('inputfile "%s" does not exist!', inputfile), end
[d, currentName, ext ] = fileparts(inputfile);
locs = file_locs(inputfile,taskdirectory,task);

% cap locations
eeglabpath = fileparts(which('eeglab'));
cap_location = fullfile(eeglabpath,'/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
if ~exist(cap_location, 'file'), error('cannot find file for 128 channel cap: %s', cap_location), end
correction_cap_location = hera('Projects/7TBrainMech/scripts/eeg/Shane/resources/ChanLocMod128to64.ced');
if ~exist(correction_cap_location, 'file'), error('cannot find file for correction 128 channel cap: %s', correction_cap_location), end

% to know how far your script is with running
fprintf('==========\n%s:\n\t Initial Preprocessing(%s,%s)\n',...
    currentName, inputfile, taskdirectory)

% where to save things
filterfolder = 'filtered';

% channel location
commonPlus = {'AFz','C1','C2','C3','C4','C5','C6','CP1','CP2','CP3','CP4',...
    'CP5','CP6','CPz','Cz','F1','F2','F3','F4','F5','F6','F7','F8','FC1',...
    'FC2','FC3','FC4','FC5','FC6','FCz','Fp1','Fp2','FT10','FT9','Fz','I1',...
    'I2','O1','O2','Oz','P1','P10','P2','P3','P4','P5','P6','P7','P8','P9',...
    'PO10','PO3','PO4','PO7','PO8','PO9','POz','Pz','T7','T8',...
    'AF8','AF7','AF4','AF3'};

if condition == 1
    xEEG = load_if_exists(locs.filter);
end

if isstruct(xEEG)
    [ALLEEG EEG CURRENTSET] = pop_newset([], xEEG, 0,...
        'setname',currentName,...
        'gui','off');
else
    %% 1. re-reference to mastoid
    EEG = pop_loadset(inputfile);

    if size(EEG.data,1) < 100
        Flag128 = 0;
        EEG = pop_reref(EEG, [65 66]); %does this "restore" the 40dB of "lost SNR" ? -was it actually lost? ...this is potentially undone by PREP
        EEG = eeg_checkset(EEG);
        % find the cap to use
    else
        %[129 130] are the mastoid externals for the 128 electrode
        EEG = pop_reref(EEG, [129 130]); %does this "restore" the 40dB of "lost SNR" ? -was it actually lost? ...this is potentially undone by PREP
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

    %% 2. filtering
    EEG = pop_eegfiltnew(EEG, lowBP, topBP, 3380, 0, [], 0);
    % filtorder = 3380 - filter order (filter length - 1). Mandatory even. performing 3381 point bandpass filtering.
    % lowBP = 0.5 Hz
    % topBP = 70 Hz

    %give a new setname and overwrite unfiltered data
    EEG = pop_editset(EEG,'setname',[currentName '_bandpass_filtered']);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    %% 3. resampling

    % Downsample the data to 150Hz using antialiasing filter
    EEG = pop_resample(EEG, 150, 0.8, 0.4); %0.8 is fc and 0.4 is dc. Default is .9 and .2. We dont know why Alethia changed them

    EEG = eeg_checkset(EEG);

    %change setname
    EEG = pop_editset(EEG,'setname',[currentName '_filtered']);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    %save filtered data
    EEG = pop_saveset( EEG, 'filename',[currentName '_filtered'], ...
        'filepath',fullfile(taskdirectory,filterfolder));
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);

    %% 4. epoching
    % skipping if already created
    if exist(locs.epoch, 'file')
        OUTEEG = [];
        fprintf('skipping; already created %s\n', locs.epoch);
        return
    end
    if ~exist(inputfile,'file')
        error('inputfile "%s" does not exist!', inputfile) 
    end


    % [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); % not sure what this is for

    EEG = pop_epoch( EEG, {'2'}, [-0.8  0.8], 'newname', [currentName '_filtered_epochs'], 'epochinfo', 'yes');

    % save epoched data set
    EEG = pop_saveset( EEG,'filename',[currentName '_filtered_epochs'], 'filepath',epochfolder);

end
end