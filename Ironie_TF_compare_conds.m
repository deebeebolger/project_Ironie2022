close all
clear all

ROInames = {'Anterior-Ventral Left', 'Anterior-Dorsal Left', 'Posterior-Ventral Left', 'Posterior-Dorsal Left',...
    'Anterior Midline', 'Posterior Midline','Anterior-Ventral Right', 'Anterior-Dorsal Right', 'Posterior-Ventral Right', 'Posterior-Dorsal Right'};
chansoi = ROInames;
bltype  = "dB";
%ROInames_order = Chanoi;
ROInames_order = {'Anterior-Ventral Left', 'Anterior Midline', 'Anterior-Ventral Right', 'Anterior-Dorsal Left', 'Posterior Midline', 'Anterior-Dorsal Right',...
    'Posterior-Ventral Left', 'Posterior-Ventral Right', 'Posterior-Dorsal Left','Posterior-Dorsal Right'};
chan2plot = [1 1 1 1 1 1 1 1 1].*5;

matfile_path = fullfile(filesep,'Users','bolger','Documents','work','Projects','Ironie','Data', ...
    'timefreq_mats', 'GA_timefreq',filesep); % Path to the folder with TF mat files.

Current_cond = {'NSI','NSL', 'Difference', 'SI', 'SL', 'Difference', 'PI', 'PL', 'Difference'};

%%

timefile = dir(fullfile(matfile_path,'*Times.mat'));
TIn = struct2cell(load(fullfile(timefile.folder,timefile.name)));
time = TIn{1,1}(1:end);
freqfile = dir(fullfile(matfile_path,'*FreqOI.mat'));
FreqIn = struct2cell(load(fullfile(freqfile.folder,freqfile.name)));
freqIn = FreqIn{1,1};

%% Set up the figure for plotting

hf1 = figure;
set(hf1,'NumberTitle', 'off', ...
    'Name', sprintf('Participant (n=%d)',22));
scr_siz = get(0,'ScreenSize') ;
pos = floor([scr_siz(3) scr_siz(4) scr_siz(3) scr_siz(4)]);
set(hf1,'Position',pos,'Color',[1 1 1]);,
row_cols = [3 3];
clim = [-0.8 0.8];
tindx = find([time>=-0.25 & time<=1.0]);


for currcond = 1:size(Current_cond,2)


    if mod(currcond, 3) == 1

        filenames    = dir(fullfile(matfile_path,strcat(Current_cond{1,currcond}, '_GAtimefreq.mat')));

        fprintf('Loading in mat-file: %s\n', filenames(1).name);
        DIn = load(fullfile(filenames(1).folder,filenames(1).name));
        CurrData = DIn.GAtimefreq;
        Diff1 = CurrData;

    elseif mod(currcond, 3) == 2

        filenames    = dir(fullfile(matfile_path,strcat(Current_cond{1,currcond}, '_GAtimefreq.mat')));

        fprintf('Loading in mat-file: %s\n', filenames(1).name);
        DIn = load(fullfile(filenames(1).folder,filenames(1).name));
        CurrData = DIn.GAtimefreq;
        Diff2 = CurrData;

    elseif mod(currcond, 3) == 0

        CurrDiff = Diff1{1, chan2plot(currcond)}(:,tindx) - Diff2{1, chan2plot(currcond)}(:,tindx);
    end


    %% Plot the result of applying the complex Morlet wavelet transform.

    if mod(currcond, 3) > 0

        fprintf('Current channel %s\n',string(ROInames{1,chan2plot(currcond)}));
        ax1(currcond) = subplot(row_cols(1),row_cols(2),currcond);
        imagesc('Parent',ax1(currcond),'XData',time(tindx),'YData',FreqIn{1,1},'CData',CurrData{1,chan2plot(currcond)}(:,tindx),'CDataMapping','scaled')
        ax1(currcond).YLim = [4 20];                    %Set Y-axis limites
        ax1(currcond).CLim = clim;                      %Set time-frequency power limits.
        ax1(currcond).XLim = [-0.25 0.8];
        ax1(currcond).Layer = 'top';
        ax1(currcond).XLabel.String = 'Time (ms)'; %Set x-axis label
        ax1(currcond).YLabel.String = 'Frequency (Hz)'; %Set y-axis label
        ax1(currcond).FontSize = 14;


        [tt,ss] = title(strcat(Current_cond{1,currcond}, ': ', ROInames{1, chan2plot(currcond)}));
        tt.FontSize = 14;
        colormap(ax1(currcond),jet);
        cb(currcond) = colorbar;                                 %Colorbar
        cb(currcond).Label.String = ['ERSP (',bltype,' )']; %Colorbar title defined
        cb(currcond).FontSize = 14;
        set(ax1(currcond),'HitTest','on','SelectionHighlight','on','UserData',{time,tindx,FreqIn{1,1},CurrData{1,chan2plot(currcond)},clim,ROInames_order{chan2plot}, bltype},'Nextplot','replace');
        set(ax1(currcond),'ButtonDownFcn',@plotsingle_tfGA)

        hold on

        f = [1 2 3 4];
        v = [0.5 4; 0.8 4; 0.8 7; 0.5 7];y
        patch('Faces', f, 'Vertices', v, 'EdgeColor', 'k', 'FaceColor', 'none', 'LineWidth', 2, 'LineStyle', '--')

    elseif mod(currcond, 3) == 0
        
        ax1(currcond) = subplot(row_cols(1),row_cols(2),currcond);
        imagesc('Parent',ax1(currcond),'XData',time(tindx),'YData',FreqIn{1,1},'CData',CurrDiff,'CDataMapping','scaled')
        ax1(currcond).YLim = [4 20];                    %Set Y-axis limites
        ax1(currcond).CLim = clim;                      %Set time-frequency power limits.
        ax1(currcond).XLim = [-0.25 0.8];
        ax1(currcond).Layer = 'top';
        ax1(currcond).XLabel.String = 'Time (ms)'; %Set x-axis label
        ax1(currcond).YLabel.String = 'Frequency (Hz)'; %Set y-axis label
        ax1(currcond).FontSize = 14;


        [tt,ss] = title(strcat(Current_cond{1,currcond}, '-Difference : ', ROInames{1, chan2plot(currcond)}));
        tt.FontSize = 14;
        colormap(ax1(currcond),jet);
        cb(currcond) = colorbar;                                 %Colorbar
        cb(currcond).Label.String = ['ERSP (',bltype,' )']; %Colorbar title defined
        cb(currcond).FontSize = 14;
        set(ax1(currcond),'HitTest','on','SelectionHighlight','on','UserData',{time,tindx,FreqIn{1,1},CurrDiff,clim,ROInames_order{chan2plot}, bltype},'Nextplot','replace');
        set(ax1(currcond),'ButtonDownFcn',@plotsingle_tfGA)

        hold on

        f = [1 2 3 4];
        v = [0.5 4; 0.8 4; 0.8 7; 0.5 7];
        patch('Faces', f, 'Vertices', v, 'EdgeColor', 'k', 'FaceColor', 'none', 'LineWidth', 2, 'LineStyle', '--')

    end




end