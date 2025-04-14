function [suberspdata,times,freqs,failcode,subID,scanDate,age,sub_ersp_data_name] = calc_ersp(epoch_path,sub_epoch_filename,sub_ersp_path,subMerge7t,channelTemplate,task)
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

% determine if subject's ersp has already been calculated
if exist(fullfile(sub_ersp_path,sub_ersp_data_name),'file')
    fprintf('computed ersp for %s\nLoading suberspdata\n',lunaid)
    subStruct = load(fullfile(sub_ersp_path,sub_ersp_data_name),'suberspdata','times','freqs');
    suberspdata = subStruct.suberspdata;
    times = subStruct.times;
    freqs = subStruct.freqs;
    EEG = pop_loadset(inputfile);
    age = EEG.age;
    failcode = NaN;
    return
else
    % skipping if subject is not in merge 7t
    if isempty(subMerge7t)
        fprintf('\nSubject %s Merge7t\n',lunaid)
        failcode = 1;
        suberspdata = [];
        times = [];
        freqs = [];
        age = [];
        return
    end
    
%     if task == "anti"
%         % determing if subject is viable
%         if subMerge7t.Viable ==0
%             fprintf('\nSubject %s not viable\n',lunaid)
%             failcode = 0;
%             suberspdata = [];
%             times = [];
%             freqs = [];
%             age=[];
%             return
%         end
%     end
    
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
        times = [];
        freqs = [];
        return
    end
    
    % skip if subject has zero epochs
    numEpochs = length(EEG.epoch);
    if numEpochs < 1
        fprintf('skipping; %s has %d epochs\n',lunaid,numEpochs)
        failcode = 2;
        suberspdata = [];
        times = [];
        freqs = [];
        return
    end
    
    % calculate subject's ersp for each channel
    for currentChan = 1:numchans
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
        % Filtered original data with bandpass filter 0.5-70 Hz -> should not look at frequencies above 70 Hz
        [suberspdata(currentChan,:,:),~,~,times,freqs]=newtimef(EEG.data(currentChan,:,:),...
             EEG.pnts,[EEG.xmin EEG.xmax]*1000, EEG.srate, [3 0.5], 'freqs', [3 70],...
            'baseline',[-700 -400],'plotersp','off','plotitc','off','plotphasesign','off'); 
    end
     % subsetting ersp data for just prep period
    idx580ms = find(times==580);
    suberspdata = suberspdata(:,:,1:idx580ms);
    times = times(1:idx580ms);
    failcode = NaN;
end