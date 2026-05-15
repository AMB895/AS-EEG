%% compile ersp and itc for each task into large 4d  matrix
% subject x freqs x times x channels
clear
% adding necessary paths
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/fieldtrip-20180926/')
addpath('/Volumes/Hera/Abby/Resources/eeglab_current/eeglab2024.2/')
addpath('/Volumes/Hera/Projects/7TBrainMech/scripts/eeg/Shane/Preprocessing_Functions/')
addpath('/Volumes/Hera/Abby/AS_EEG/')
addpath('/Volumes/Hera/Abby/Resources/slanCM/')

%% Setting up directories
% merged 7t  eeg info and AS scores
merged7t_eeg = readtable('/Volumes/Hera/Abby/AS_EEG/Results/merged7t_ASscore_table.csv');
tasks = ["corAS","VGS","errcorAS"];
eventtype = 'prep';
basenorm = 'avgtrl';
% basenorm = 'singletrl';
maindir = hera('Abby/AS_EEG/ERSPdata/');

% Rerun?
RERUN=0;
% plot?
PLOT=1;
SAVE_FROW = 1;


% loop of tasks
for t = 1:length(tasks)
    currtask = char(tasks(t));
    taskdir_subs = [maindir currtask '_data/' eventtype '/SubjectData/' basenorm '/'];
    subfiles = dir(fullfile(taskdir_subs,'*.mat'));
    % define save path
    taskdir_ersp = sprintf('/Volumes/Hera/Abby/AS_EEG/ERSPdata/%s_data/%s',currtask,eventtype);
    % file name
    ersp_name = sprintf('%s_%s_%s_ersp.mat',currtask,eventtype,basenorm);
    itc_name = sprintf('%s_%s_%s_itc.mat',currtask,eventtype,basenorm);
    idmat_name = sprintf('%s_%s_%s_ids.mat',currtask,eventtype,basenorm);
    % waitbar
    w=waitbar(0,"Starting","Name",sprintf("All ERSP %s %s",currtask,eventtype));
    
    % skip if already saved
    if exist(fullfile(taskdir_ersp,idmat_name),'file') && RERUN==0
        fprintf('computed ERSP and time/freq data for %s %s\n',currtask,eventtype)
        close(w)
        if PLOT 
            MAP = slanCM('viridis',256);
            % load total ersp
            load(fullfile(taskdir_ersp,ersp_name))
            load(fullfile(taskdir_ersp,itc_name))
            load(fullfile(taskdir_ersp,idmat_name))
            % remove non viable
            erspdata = squeeze(mean(allersp(idmat(:,4)==1,:,:,[4 5 6 7 37 38 39 40 41]),[4]));
            itcdata = squeeze(mean(allitc(idmat(:,4)==1,:,:,[4 5 6 7 37 38 39 40 41]),[4]));
            if SAVE_FROW
                avg_ersp = squeeze(mean(erspdata,1));
                avg_itc = squeeze(mean(itcdata,1));
                save(sprintf("/Volumes/Hera/Abby/AS_EEG/ERSPdata/%s_data/prep/%s_prep_avgtrl_frow_ersp.mat",currtask,currtask),"avg_ersp","times","freqs")
                save(sprintf("/Volumes/Hera/Abby/AS_EEG/ERSPdata/%s_data/prep/%s_prep_avgtrl_frow_itc.mat",currtask,currtask),"avg_itc","times","freqs")
            end
            idmat_viable = idmat(idmat(:,4)==1,:);
            erspdata_adults = squeeze(mean(erspdata(idmat_viable(:,3)>=18,:,:),1));
            itcdata_adults = squeeze(mean(itcdata(idmat_viable(:,3)>=18,:,:),1));
            erspdata_adolescents = squeeze(mean(erspdata(idmat_viable(:,3)<18,:,:),1));
            itcdata_adolescents = squeeze(mean(itcdata(idmat_viable(:,3)<18,:,:),1));
            
            % plot average f-row itc
%             plot_tfce_clusters(times,freqs,squeeze(mean(erspdata,1)),...
%                 'title',[currtask ' ERSP ' basenorm],'subtitle','F-row electrodes',...
%                 'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','dB','addlines',1,'map',MAP)
%             % plot average f-row itc just adults
%             plot_tfce_clusters(times,freqs,erspdata_adults,...
%                 'title',[currtask ' ERSP adults ' basenorm],'subtitle','F-row electrodes',...
%                 'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','dB','addlines',1,'map',MAP)
%             % plot average f-row itc just adults
%             plot_tfce_clusters(times,freqs,erspdata_adolescents,...
%                 'title',[currtask ' ERSP adolescents ' basenorm],'subtitle','F-row electrodes',...
%                 'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','dB','addlines',1,'map',MAP)
%             
%             % plot average f-row itc
%             plot_tfce_clusters(times,freqs,squeeze(mean(itcdata,1)),...
%                 'title',[currtask ' ITC'],'subtitle','F-row electrodes',...
%                 'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','Coherence','addlines',1,'map',MAP)
%             % plot average f-row itc just adults
%             plot_tfce_clusters(times,freqs,itcdata_adults,...
%                 'title',[currtask ' ITC adults'],'subtitle','F-row electrodes',...
%                 'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','Coherence','addlines',1,'map',MAP)
%             % plot average f-row itc just adolescents
%             plot_tfce_clusters(times,freqs,itcdata_adolescents,...
%                 'title',[currtask ' ITC adolescents'],'subtitle','F-row electrodes',...
%                 'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','Coherence','addlines',1,'map',MAP)
        end
    else
        % set up large matricies
        allersp = zeros(length(subfiles),150,76,64);
        idmat = zeros(length(subfiles),4);
        allitc = zeros(length(subfiles),150,76,64);
        
        % loop over subject files
        for currsub = 1:length(subfiles)
            % get lunaid and eeg date
            filename = subfiles(currsub).name;
            splitname = split(filename,"_");
            lunaid = splitname{1};
            eeg_date = splitname{2};
            % update waitbar
            waitbar(currsub/length(subfiles),w,sprintf("%s %s",lunaid,eeg_date))
            % defining subject's merged7t table
            subIdx_7t = find(merged7t_eeg.lunaid==str2double(lunaid) & merged7t_eeg.eeg_date==str2double(eeg_date));
            subTable_7t = merged7t_eeg(subIdx_7t,:);
            eeg_age = subTable_7t.eeg_age;
            % add to id, ersp, and itc matricies
            idmat(currsub,:) = [str2double(lunaid), str2double(eeg_date), eeg_age, subTable_7t.Viable];
            % load subject's data
            load(fullfile(taskdir_subs,filename),'sub_ersp','sub_itc','times','freqs');
            allersp(currsub,:,:,:) = sub_ersp;
            allitc(currsub,:,:,:) = sub_itc;
        end
        % save idmat, allersp, and allitc
        save(fullfile(taskdir_ersp,ersp_name),'allersp','times','freqs')
        save(fullfile(taskdir_ersp,itc_name),'allitc','times','freqs')
        save(fullfile(taskdir_ersp,idmat_name),'idmat')
        close(w)
        
        if PLOT
            MAP = slanCM('viridis',256);
            % load total ersp
            load(fullfile(taskdir_ersp,ersp_name))
            load(fullfile(taskdir_ersp,itc_name))
            load(fullfile(taskdir_ersp,idmat_name))
            % remove non viable
            erspdata = squeeze(mean(allersp(idmat(:,4)==1,:,:,[4 5 6 7 37 38 39 40 41]),[4]));
            itcdata = squeeze(mean(allitc(idmat(:,4)==1,:,:,[4 5 6 7 37 38 39 40 41]),[4]));
            idmat_viable = idmat(idmat(:,4)==1,:);
            erspdata_adults = squeeze(mean(erspdata(idmat_viable(:,3)>=18,:,:),1));
            itcdata_adults = squeeze(mean(itcdata(idmat_viable(:,3)>=18,:,:),1));
            erspdata_adolescents = squeeze(mean(erspdata(idmat_viable(:,3)<18,:,:),1));
            itcdata_adolescents = squeeze(mean(itcdata(idmat_viable(:,3)<18,:,:),1));
            
            % plot average f-row itc
            plot_tfce_clusters(times,freqs,squeeze(mean(erspdata,1)),...
                'title',[currtask ' ERSP ' basenorm],'subtitle','F-row electrodes',...
                'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','dB','addlines',1,'map',MAP)
            % plot average f-row itc just adults
            plot_tfce_clusters(times,freqs,erspdata_adults,...
                'title',[currtask ' ERSP adults ' basenorm],'subtitle','F-row electrodes',...
                'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','dB','addlines',1,'map',MAP)
            % plot average f-row itc just adults
            plot_tfce_clusters(times,freqs,erspdata_adolescents,...
                'title',[currtask ' ERSP adolescents ' basenorm],'subtitle','F-row electrodes',...
                'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','dB','addlines',1,'map',MAP)
            
            % plot average f-row itc
            plot_tfce_clusters(times,freqs,squeeze(mean(itcdata,1)),...
                'title',[currtask ' ITC'],'subtitle','F-row electrodes',...
                'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','Coherence','addlines',1,'map',MAP)
            % plot average f-row itc just adults
            plot_tfce_clusters(times,freqs,itcdata_adults,...
                'title',[currtask ' ITC adults'],'subtitle','F-row electrodes',...
                'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','Coherence','addlines',1,'map',MAP)
            % plot average f-row itc just adolescents
            plot_tfce_clusters(times,freqs,itcdata_adolescents,...
                'title',[currtask ' ITC adolescents'],'subtitle','F-row electrodes',...
                'xlab','Preparatory Period (ms)','ylab','Frequency (Hz)','clab','Coherence','addlines',1,'map',MAP)
        end
        % clear variables for next task
        allersp = []; allitc = []; idmat = [];
    end
end