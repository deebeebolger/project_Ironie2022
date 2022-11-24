%%*************************************************************************************************************
% Date: 19-09-2022            Programmed by: d. Bolger
% Script to take load in time-frequency data from the *timefreq.mat files
% and to average TF data over defined frequency bands and time-windows.
% The data needs to be saved in an excel file.
% This is carried out for each participant.

function [datastruct, tf_powerdB, time] = TF2Excel()

    allconds = {'PL', 'PI', 'NSL', 'NSI', 'SL', 'SI'};
    condfull = {'Control_Literal', 'Control_Ironic', 'NonSarcastic_Literal', 'NonSarcastic_Ironic', 'Sarcastic_Literal', 'Sarcastic_Ironic'};
    freqbands_nom = {'theta', 'alpha1', 'alpha2', 'beta1','beta2', 'gamma1', 'gamma2'};
    freqbands_hz = {[4 7], [8 10], [11 13], [14 20], [21 30], [30 45], [60 80]};
    
    ROIindx = [1, 7; 8, 14; 15 20; 21, 26; 27, 30; 31, 34; 35, 41; 42, 48; 49, 54; 55, 60];
    ROInames = {'AV_left', 'AD_left', 'PV_left', 'PD_left', 'AMid', 'MPos','AV_right', 'AD_right', 'PV_right', 'PD_right'};
    chanois = {'AF7', 'F7', 'FT7', 'F5', 'FC5', 'C5', 'T7','AF3', 'F3', 'FC3', 'F1', 'FC1', 'C1', 'C3','TP7', 'CP5', 'P5', 'P7', 'P9', 'PO7',...
    'CP3', 'CP1', 'P3', 'PO3', 'P1', 'O1','AFz', 'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'POz','Oz','AF8', 'F8', 'FT8', 'F6', 'FC6', 'C6', 'T8','AF4', 'F4', 'FC4', 'F2', 'FC2', 'C2', 'C4',...
    'TP8', 'CP6', 'P6', 'P8', 'P10', 'PO8','CP4', 'CP2', 'P4', 'PO4', 'P2', 'O2'};

    %% Need to find the indices of the channels in the 64-channel context.
    chaninfo = load("/Users/bolger/Documents/work/Projects/Ironie/Data/Chaninfo.mat");
    chans_all = {chaninfo.chaninfo.labels};

    ROIs64_indx = cell(size(ROInames, 2), 1);

    for rcnt = 1:size(ROInames,2)
        ROIs64_indx{rcnt,1} = find(ismember(chans_all, chanois(ROIindx(rcnt,1):ROIindx(rcnt,2))));
    
    end
    
    %% Present a dialogue box so that the user can entre the temporal window
    % limits (in seconds) over which the data will be averaged.
    prompt = {'Enter the min and max in seconds of each time window','Name the time windows'};
    dlgtitle = 'Define time windows';
    dims     = [10 25];
    tans     = inputdlg(prompt, dlgtitle, dims);
    time_winds = string(tans{1,1});
    twind_name = cellstr(tans{2,1});
    
    Twinds = cellfun(@split, time_winds, 'UniformOutput',false);   % stored as a cell array of strings.
    
    datastruct = [];
    
    for icond = 1:length(allconds)
    
        % Extract the current occupation-type and speech-type.
        currcond = allconds{1,icond};
        C = split(condfull{1,icond},'_');
        fprintf('The current occupation is %s and current speech-type is %s\n', C{1,1}, C{2,1})
        occupation_type = C{1,1};   % Control, NonSarcastic, Sarcastic
        speech_type    = C{2,1};   % Literal, Ironic
    
        % Define the path to current condition-specific, participant-level *.mat file.
       
        matfile_path = fullfile(filesep,'Users','bolger','Documents','work','Projects','Ironie','Data','timefreq_mats', currcond,filesep); % Path to the folder with TF mat files.
        filenames    = dir(fullfile(matfile_path, '*timefreq.mat'));
        allnames     = {filenames.name};
    
        % Initialize variables.
        DataIn = cell(1,numel(filenames));
        erspIn = cell(1, numel(filenames));
        freqIn = cell(1, numel(filenames));
    
        % Only need to load in the time vector once (at the start).
        if icond ==1
            timefile = dir(fullfile(matfile_path,'Times.mat'));         % Time vector is in msec so need to convert to seconds.
            TIn = struct2cell(load(fullfile(timefile.folder,timefile.name), 'Tnew'));
            time = TIn{1,1}(1:end);
    
        else
            fprintf('Time vector is already loaded!')
        end
    
    
        for i = 1:numel(filenames)                % Load in the subject-level TF mat files individually.
    
            fprintf('Loading in mat-file: %s\n', filenames(i).name)
            DataIn{1,i} = load(fullfile(filenames(i).folder,filenames(i).name));
            currdata    = DataIn{1,i}.timefreq_results;
            erspIn{1,i} = currdata.ersp;
            freqIn{1,i} = currdata.freqs;
            %chans = DataIn{1,1}.timefreq_results.chansoi{1,1};

            erspData = cell(size(freqbands_hz,2), numel(filenames));
            ROImean  = cell(size(freqbands_hz,2), numel(filenames));
    
            for freqcnt = 1:size(freqbands_hz,2)
    
                % Find the indices of data for the current frequency band of
                % interest
    
                fcurr_indx = dsearchn(freqIn{1,i}', freqbands_hz{1,freqcnt}');
                erspData{freqcnt,i} = getcurr_window(erspIn{1,i}, time, Twinds, fcurr_indx);  % Call of function to extract TF data for the current subject.

                % Here function to average channel ROIs
                ROImean{freqcnt,i} = ROI_calc(erspData{freqcnt,i}, chanois, ROIs64_indx, ROInames);
    
                datastruct.(occupation_type).(speech_type).(strcat('S',allnames{1,i}(1:end-22))).(freqbands_nom{1,freqcnt}) = ROImean{freqcnt,i}';
            end
    
        end
    
    end   % end of icond loop
    assignin('base','datastruct',datastruct)
    assignin('base','allnames', allnames)
    assignin('base', 'time', time)
    
    % Order of data: Occupation-type, Speech-type (Ironic, literal),
    % participant-title, channel, frequency-band, time-wind1, time-wind2,
    % time-wind3
    % Create the Subject Column
    datatable = table();
    A = squeeze(split(condfull,'_'));
    assignin('base','A',A)
    occup_all = unique({A{:,1}});
    speech_all = unique({A{:,2}});
    savexl_path = fullfile(filesep,'Users','bolger','Documents','work','Projects','Ironie','Data',filesep);
    xlname = 'Ironie_TFdata_stats_ROI_latewindow_v4.xlsx';
    
    
    for sujcntr = 1:size(allnames,2)
    
        currsuj_cols = size(ROInames,2) * size(Twinds,1) * size(freqbands_nom,2) * 2 * 3;
    
        currsuj_nom = cellstr(allnames{1,sujcntr}(1:end-22));
        Subject = repmat(currsuj_nom, currsuj_cols, 1);
        assignin('base', "Subject", Subject)
    
        tf_powerdB = [];
        for ioccup = 1:size(occup_all,2)
            for ispeech = 1:size(speech_all,2)
                for ifreq = 1:size(freqbands_hz,2)
                    twind_data = datastruct.(occup_all{1,ioccup}).(speech_all{1,ispeech}).(strcat('S',allnames{1,sujcntr}(1:end-22))).(freqbands_nom{1,ifreq});
                     tf_powerdB = cat(1, tf_powerdB,reshape(twind_data, [size(twind_data,2)*size(twind_data,1),1]));
                end
            end
        end
        assignin('base', 'tf_powerdB', tf_powerdB)
        

        chanrep = currsuj_cols/size(ROInames,2);
        Channels = repmat(ROInames', chanrep,1);
        assignin('base', 'Channels', Channels)

        TFw1 = reshape(repmat(twind_name', size(ROInames,2), 1), [size(ROInames,2)*size(twind_name,1),1]);
        twindrep = currsuj_cols/length(TFw1);
        Time_window = reshape(repmat(TFw1, twindrep, 1),[twindrep*size(TFw1,1),1]);
        assignin('base', 'Time_window', Time_window)

        FB1 = reshape(repmat(freqbands_nom,size(TFw1,1),1),[size(TFw1,1)*size(freqbands_nom,2),1]);
        freqrep = currsuj_cols/size(FB1,1);
        FrequencyBand = reshape(repmat(FB1, freqrep,1),[freqrep*size(FB1,1),1]);
        assignin('base', 'FrequencyBand', FrequencyBand)

        Spch1 = reshape(repmat(speech_all, size(FB1,1),1),[size(FB1,1)*size(speech_all,2),1]);
        speechrep = currsuj_cols/size(Spch1,1);
        Annonce_type = reshape(repmat(Spch1, speechrep,1), [speechrep*size(Spch1,1),1]);
        assignin('base', 'Annonce_type', Annonce_type)

        Occup1 = reshape(repmat(occup_all, size(Spch1,1),1),[size(Spch1,1)*size(occup_all,2),1]);
        occuprep = currsuj_cols/size(Occup1,1);
        Occupation_type = reshape(repmat(Occup1, occuprep,1), [occuprep*size(Occup1,1),1]);
        assignin('base', 'Occupation_type', Occupation_type)

       
        currtable = table(Subject, Occupation_type, Annonce_type, FrequencyBand, Time_window, Channels, tf_powerdB);
        assignin('base', 'currtable', currtable)
        writetable(currtable, fullfile(savexl_path,xlname), 'Sheet',strcat('S',allnames{1,sujcntr}(1:end-22)), 'WriteVariableNames',true);

    end

end   % end of main function 



function ersp_wind = getcurr_window(tfIn, T, twinds, findx)
% Function to extract the time-window data from the ersp data for the
% current participants.
% Input:
% - tfIn: input time-frequency data for current frequency band and
% participant (1 X no. channels cell array)
% - T: time vector (1 X time-points)
% - twindw: time window limits (no. time windows X 1).
% - findx : indices of upper and lower limits of current frequency band.
% Output: ersp_wind = cell array of tf data for current frequency band
% (dimensions: number of time-window X number of channels).
%*******************************************************************************


ersp_wind = cell(size(twinds,1), size(tfIn,2));
for tcnt = 1:size(twinds,1)

    Tcurr      = twinds{tcnt, 1};
    Tlims      = [str2double(Tcurr{1,1}), str2double(Tcurr{2,1})];
    Tlims_indx = dsearchn(T'*1000, Tlims');

    for ichan = 1:size(tfIn,2)

        currTF = tfIn{1,ichan};
        fprintf('size of currTF is %d\n', size(currTF))
        ersp_wind{tcnt, ichan} = mean(mean(currTF(findx(1):findx(2),Tlims_indx(1):Tlims_indx(2)),2),1);  % Average over frequency-band & time-window.

    end
end
ersp_windmat = cell2mat(ersp_wind)

end % end of function twind_data()

%%----------- Call of subfunction to calculation mean over ROIs ------------------

function tfmean = ROI_calc(tfdataIn, Chans, roi_indx, roi_names)

tdata = cell2mat(tfdataIn);
tfmean = zeros(size(tdata,1), size(roi_indx,1));
for icnt = 1:size(roi_indx,1)
    tfmean(:,icnt) = mean(tdata(:,roi_indx{icnt,1}),2);
end


end
    


