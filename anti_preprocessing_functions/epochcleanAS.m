function [revisar] = epochcleanAS(inputfile,epoch_folder,kept_epoch_folder,task,taskdirectory,trialtype)
revisar = {};
% what file are we using

% locs = file_locs(inputfile,epoch_folder,task);
% locs = file_locs(inputfile,taskdirectory,task);
% if exist(locs.epochClean, 'file')
%     OUTEEG = [];
%     fprintf('skipping; already created %s\n', locs.epochClean);
%     %OUTEEG = pop_loadset(locs.ICAwholeClean);
%     return
% end

if ~exist(inputfile,'file')
    error('inputfile "%s" does not exist!', inputfile) 

end

[d, currentName, ext ] = fileparts(inputfile);

if trialtype == 1
    epoch_name = [currentName '_epochs_cor'];
    epochmarked_name = [currentName '_epochs_marked_cor'];
    epochkept_name = [currentName '_epochs_kept_cor'];
    if exist(epochkept_name,'file')
        OUTEEG = [];
        fprintf('skipping; already created %s\n',epochkept_name)
        return
    end
elseif trialtype == 0
    epoch_name = [currentName '_epochs_incor'];
    epochmarked_name = [currentName '_epochs_marked_incor'];
    epochkept_name = [currentName '_epochs_kept_incor'];
    if exist(epochkept_name,'file')
        OUTEEG = [];
        fprintf('skipping; already created %s\n',epochkept_name)
        return
    end
elseif trialtype == 2
    epoch_name = [currentName '_epochs_errcor'];
    epochmarked_name = [currentName '_epochs_marked_errcor'];
    epochkept_name = [currentName '_epochs_kept_errcor'];
    if exist(epochkept_name,'file')
        OUTEEG = [];
        fprintf('skipping; already created %s\n',epochkept_name)
        return
    end
end

% where to find eeglab stuff
eeglabpath = fileparts(which('eeglab'));
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

EEG = pop_loadset(inputfile);
if EEG.nbchan~=64
    revisar = inputfile;
else
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


%% epoching: -0.500 to 2 seconds
%extract epochs
if trialtype == 1
     EEG = pop_epoch( EEG, {'2_cor'}, [-1  2], 'newname',  epoch_name, 'epochinfo', 'yes');
elseif trialtype == 0
    EEG = pop_epoch( EEG, {'2_incor'}, [-1  2], 'newname', epoch_name, 'epochinfo', 'yes');
elseif trialtype == 2
    EEG = pop_epoch( EEG, {'2_errcor'}, [-1  2], 'newname', epoch_name, 'epochinfo', 'yes');
end

% %remove baseline
% EEG = pop_rmbase( EEG, [-400    0]);

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'overwrite','on','gui','off');

%save epoched eegsets
EEG = pop_saveset( EEG,'filename',[currentName epoch_name], ...
    'filepath',epoch_folder);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% ~10% should be rejected.

%Two options for epoch rejection: 2nd option is used.

%1. kurtosis
%kurtosis (not recommended): default 5 for maximum threshold limit. Try that, check how many
%epochs removed, otherwise higher to 8-10.
% EEG = pop_autorej(EEG, 'nogui','on','eegplot','off');

%2.Use improbability and thresholding
%Apply amplitude threshold of -500 to 500 uV to remove big
% artifacts(don't capture eye blinks)
EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan],-500,500,0,EEG.xmax,0,1);


%apply improbability test with 6SD for single channels and 2SD for all channels,
EEG = pop_jointprob(EEG,1,[1:EEG.nbchan],6,2,0,0,0,[],0);

%save marked epochs (to check later if you agree with the removed opochs)
EEG = pop_editset(EEG,'setname',[currentName epochmarked_name]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_saveset( EEG,'filename',[currentName epochmarked_name], ...
    'filepath',epoch_folder);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


%reject epochs
EEG = pop_rejepoch(EEG, find(EEG.reject.rejjp), 0);

%save epochs rejected EEG data
%this saves the non-rejected epochs
EEG = pop_editset(EEG, 'setname', epochkept_name);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_saveset(EEG, 'filename', epochkept_name, 'filepath',kept_epoch_folder);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);


