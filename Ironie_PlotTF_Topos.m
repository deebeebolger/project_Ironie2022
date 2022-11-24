%%
% Programmed By: D. Bolger                 Date: 02-11-2022
% Plot the topographies of time-frequency results over a specific time
% window.
%****************

currcond = {'Literal','Ironic', 'Difference'};
toi = [500 800];
foi = [14 20];
chaninfo = load("/Users/bolger/Documents/work/Projects/Ironie/Data/Chaninfo.mat");

hf1 = figure;
set(hf1,'NumberTitle', 'off', ...
    'Name', sprintf('Participant (n=22)'));
scr_siz = get(0,'ScreenSize') ;
pos = floor([scr_siz(3) scr_siz(4) scr_siz(3) scr_siz(4)]);
set(hf1,'Position',pos,'Color',[1 1 1]);
rows = 1;
cols = 3;

%% Define the channels to plot
ROIcurr =  cell(rows, 1);
ROIcurr{1,1}  =  '';
ROIcurr{2,1}  =  '';
ROIcurr{3,1}  =  '';

ROIindx  = [1, 7; 8, 14; 15 20; 21, 26; 27, 30; 31, 34; 35, 41; 42, 48; 49, 54; 55, 60];
ROInames = {'AV_left', 'AD_left', 'PV_left', 'PD_left', 'AMid', 'MPos','AV_right', 'AD_right', 'PV_right', 'PD_right'};
chanois  = {'AF7', 'F7', 'FT7', 'F5', 'FC5', 'C5', 'T7','AF3', 'F3', 'FC3', 'F1', 'FC1', 'C1', 'C3','TP7', 'CP5', 'P5', 'P7', 'P9', 'PO7',...
    'CP3', 'CP1', 'P3', 'PO3', 'P1', 'O1','AFz', 'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'POz','Oz','AF8', 'F8', 'FT8', 'F6', 'FC6', 'C6', 'T8','AF4', 'F4', 'FC4', 'F2', 'FC2', 'C2', 'C4',...
    'TP8', 'CP6', 'P6', 'P8', 'P10', 'PO8','CP4', 'CP2', 'P4', 'PO4', 'P2', 'O2'};
chans_all = {chaninfo.chaninfo.labels};
ROIs_indx = cell(size(ROIcurr,1), 1);


for rcnt1 = 1:size(ROIcurr,1)

    curr_roi = string(ROIcurr{rcnt1, 1});
    if ~strcmp(curr_roi, '')
        roi_inc   = [];
        for rcnt = 1:size(curr_roi,2)
            rindx         = find(strcmp(ROInames, curr_roi{1,rcnt}));
            roi_indx_curr = ROIindx(rindx,:);
            R             = find(ismember(chans_all, chanois(roi_indx_curr(1):roi_indx_curr(2))));
            roi_inc       = cat(2, roi_inc, R);
        end
        ROIs_indx{rcnt1, 1} = roi_inc;
    else
        ROIs_indx{rcnt1, 1} = [];
    end
end


R2plot_indx = sort(reshape(repmat([1:rows]', 1, cols),[], rows*cols));
T_GA       = cell(numel(currcond),1);

for condcount = 1:size(currcond,2)

    if strcmp(currcond{1,condcount}, 'Difference')
        T_GA{condcount,1}   = T_GA{condcount-2,1} - T_GA{condcount-1,1};
    else
       

        matfile_path = fullfile(filesep,'Users','bolger','Documents','work','Projects','Ironie','Data', ...
            'timefreq_mats', currcond{1, condcount},filesep); % Path to the folder with TF mat files.
        % matfile_path = fullfile(filesep,'Users','bolger','Documents','work','Projects','Projet-L2-VREEG','TFmats_Test',filesep); % Path to the folder with TF mat files.
        filenames    = dir(fullfile(matfile_path,'*timefreq.mat'));


        DataIn = cell(1,numel(filenames));
        powIn = cell(1, numel(filenames));
        baselineIn = cell(1, numel(filenames));
        freqIn     = cell(1, numel(filenames));
        baseline_all = cell(1, numel(filenames));

        blup_lim   = 0;
        timefile   = dir(fullfile(matfile_path,'*Times.mat'));
        TIn        = struct2cell(load(fullfile(timefile.folder,timefile.name)));
        time       = TIn{1,1}(1:end);
        time_str   = num2str(time(end))
        if length(time_str)<3
            time = time*1000;
        end
        postim_idx = find(time<= 1000);
        time_int   = time(postim_idx);

        for i = 1:numel(filenames)

            fprintf('Loading in mat-file: %s\n', filenames(i).name)
            DataIn{1,i}     = load(fullfile(filenames(i).folder,filenames(i).name));
            currdata        = DataIn{1,i}.timefreq_results;
            powIn{1,i}      = currdata.meanTF;                % Choose the non-baseline corrected data.
            baselineIn{1,i} = currdata.baseline;              % Baseline data.
            freqIn{1,i}     = currdata.freqs;                 % The frequency vector
            erspIn{1,i}     = currdata.ersp;
        end

        lowlim = dsearchn(time', toi(1));   % Find the indices of the time-window of interest.
        uplim  = dsearchn(time', toi(2));

        lowlim_f = dsearchn(freqIn{1,1}', foi(1));
        uplim_f  = dsearchn(freqIn{1,1}', foi(2));

        TFmean_curr = cell(size(erspIn,2), 1);


        for sujcnt = 1:size(erspIn, 2)

            E         = erspIn{1, sujcnt};
            Etime_avg = cellfun(@(x) mean(x(:, lowlim:uplim),2), E, 'UniformOutput',false);
            Efreq_avg = cellfun(@(x1) mean(x1(lowlim_f:uplim_f, 1),1), Etime_avg, 'UniformOutput',false);
            TFmean_curr{sujcnt,1} = cell2mat(Efreq_avg);
        end
        T_all              = cell2mat(TFmean_curr);
        T_GA{condcount,1}  = mean(T_all,1);
    end

    %% Plot the current topomap
    ROI2plot = ROIs_indx{R2plot_indx(condcount),1};
    ax1(condcount) = subplot(rows,cols,condcount);
    topoplot(T_GA{condcount,1},chaninfo.chaninfo,'style','both','electrodes','on', 'emarker2', {ROI2plot, 'o','k'}, 'whitebk', 'on', ...
        'shading','interp', 'headrad', 'rim', 'maplimits', [-0.5 0.5])
    colorbar(ax1(condcount))


end


