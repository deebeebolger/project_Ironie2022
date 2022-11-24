%% Script to carry out time-frequency analysis of segmented EEG data.
% The time-frequency analysis is carried out using complex Morlet wavelets.
% The spectral bandwidth is defined by defining the fwhm and the empirical
% fwhm is estimated for each frequency of interest and expressed in time
% and frequency domain. 
% The resulting time-frequency map is converted to ERSP by baseline
% correction. 
% Baseline correction is carried out by calling the function
% CREX_TF_baseline(). This function gives the possibility of carry
% the following baseline corrections:
% 1. decibel conversion
% 2. Percentage change
% 3. Express post-stimulus activity as z-score in relation to baseline
% interval.
% This script carries out the time-frequency decomposition for a sin
% subject.
% This script carries out baseline decomposition for a single subject. The
% time-frequency matrix (ERSP) is saved as a matfile using the current
% filepath defined for the current dataset.
% The script expects a segmented EEGLAB *.set file.
% The ERSP is plotted for all electrodes of interest. 
% Programmed by: D. Bolger                 Date: 28-09-2021
%***************************************************************************
close all
clear all
%% Open EEGLAB 

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;               
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_loadset();
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off'); % current set = xx;
EEG = eeg_checkset( EEG );
EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',EEG.filepath);  % Saves a copy of the current resampled dataset to the current directory
eeglab redraw

%% Define the parameters.
for scounter = 1:size(ALLEEG,2)

    
    DIn = ALLEEG(scounter).data;
    T = ALLEEG(scounter).times./1000;                        % Change times from ms to seconds.
    %chansoi = {'FC3', 'FCz','FC4','C3','Cz','C4','CP3','CPz','CP4' };
    allchans = {ALLEEG(scounter).chanlocs.labels};
    chansoi  =  {allchans{1:64}};
    chanindx = find(ismember(allchans, chansoi));
    eindx = chanindx;                               % Electrode indices.
    chans = chansoi;
    
    %% Set wavelet parameters.
    
    wavet = -2:1/ALLEEG(scounter).srate:2;                   % The wavelet time vector.
    midp = dsearchn(wavet',0);                  % Find the midpoint of wavelet length.
    datalen = size(DIn,2);                      % data length (length of trial).
    wavlen = size(wavet,2);                     % wavelet length.
    nConv = wavlen+datalen-1;                   % convolution length.
    wavlen_half = floor(wavlen/2)+1;                         % half wavelet length.
    freqs = ALLEEG(scounter).srate*(0:(datalen/2))/datalen;  % Create the frequency vector.
    foi = freqs(freqs>=4 & freqs<=30);          % Define the frequency band of interest.
    
    % Define the fwhm as a function of frequency of interest.
    fwhm = linspace(0.4, 0.1, numel(foi));
    
    %% Calculate the complex Morlet wavelet for each center frequency of interest.
    
    % initialize time intervals of interest in seconds.
    tindx = [T>=-0.25 & T<=2.0];
    blindx = [T>=-0.25 & T<=0];
    Tnew = T(tindx);
    tpostindx = [Tnew>=0 & Tnew<=2.0];
    
    % Initialize variables
    tf = zeros(length(foi),length(T),size(DIn,3));
    bl = zeros(length(foi),length(find(blindx)),size(DIn,3));
    meanTF = cell(1,numel(eindx));
    meanBL = cell(1, numel(eindx));
    meanTFBL = cell(1, numel(eindx));
    
    empfwhm = zeros(length(foi),size(DIn,3),2);
    empfwhmt = zeros(length(foi),size(DIn,3),2);
    tmr = zeros(1,size(DIn,3));
    
    wb = waitbar(0,'Ready...','Name','Calculating time-frequency...',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(wb,'canceling',0);
    
    tStart = tic; 
    for ecnt = 1:length(eindx)
        tic;
        for trcnt = 1:size(DIn,3)
            
            for fcnt = 1:length(foi)
                
                Dcurr = DIn(eindx(ecnt),:,trcnt);
                DIn_fft = fft(Dcurr,nConv, 2);                                          % fft of current signal (single trial).
                
                % Define the gaussian in the time domain.
                gwin = exp( (-4*log(2)*wavet.^2) ./ fwhm(fcnt)^2 );                     % Calculate the gaussian using the current fwhm value (full-width at half maximum).
                empfwhm(fcnt,trcnt,1) = wavet(midp-1+dsearchn(gwin(midp:end)',.5));     % Calculate the empirical fwhm (time domain)
                empfwhm(fcnt,trcnt,2) = 1/empfwhm(fcnt,1);                              % Calculate the empirical fwhm (spectral domain)
                
                waveX = fft( exp(2*1i*pi*foi(fcnt)*wavet).*gwin,nConv );                % Calculate the complex Morlet wavelet for fois, given FWHM
                waveX_norm = waveX./max(waveX);                                         % Normalize
                
                DInConv = ifft(waveX_norm.*DIn_fft);                                    % Convolve (find ifft of product signal and wavelet).
                tf_pow = abs(DInConv).^2;                                               % Calculate the power
                tf(fcnt,:,trcnt) = tf_pow(wavlen_half:end-wavlen_half+1);               % Trim and reshape
                bl(fcnt,:,trcnt) = tf(fcnt,blindx,trcnt);                               % Extract the baseline.
                
            end
            
        end
        mean_tfpow = squeeze(mean(tf,3));
        mean_tfpow = mean_tfpow(:,tindx);
        meanTF{1,ecnt} = mean_tfpow;
        
        mean_blpow = squeeze(mean(bl,3));
        meanBL{1,ecnt} = mean_blpow;
        
        % Call of function to carry out baseline correction. Using the mean of
        % baseline.
        [tfpow_blc, bltype] = CREX_TF_baseline(mean_tfpow,mean_blpow,'dbels');
        meanTFBL{1,ecnt} = tfpow_blc;
        
        % Update waitbar and message
        waitbar(ecnt/length(eindx),wb,sprintf('%s',chans{1,ecnt}))
        tmr(trcnt) = toc;
    end
    
    ttotal = sum(tmr);
    disp(ttotal)
    delete(wb)
    
    %% Save the current time-frequency matrix in a structure as a mat-file.
    
    s = strfind(ALLEEG(scounter).setname,'-');
    matname = [ALLEEG(scounter).setname(1:s(1)-1),'-allchans-timefreq.mat'];
    
    timefreq_results = [];
    timefreq_results.subject = ALLEEG(scounter).setname(1:s(1)-1);
    timefreq_results.freqs = foi;
    timefreq_results.meanTF = meanTF;
    timefreq_results.ersp = meanTFBL;
    timefreq_results.baseline = meanBL;
    timefreq_results.bltype = bltype;
    timefreq_results.chansoi = {chansoi};
    
    save(fullfile(ALLEEG(scounter).filepath,matname),'timefreq_results'); % Save the time-frequency data to a mat-file in current subject folder.
    
    %% Carry out continuous wavelet transform (to get the COI)
    
    [wt,f,coi] = cwt(squeeze(DIn(48,:,8)),'amor',ALLEEG(scounter).srate);    % Based on Cz
    hf = figure;
    AX = axes('parent',hf);
    imagesc('Parent',AX,'XData',T,'YData',f,'CData',abs(wt),'CDataMapping','scaled');
    AX.Layer = 'top';
    AX.YDir = 'normal';
    AX.YScale = 'log';
    AX.XLim = [T(1) T(end)];
    hold on
    plot(AX,T,coi,'w--','linewidth',2);
    
    %% Plot the result of applying the complex Morlet wavelet transform.
    
    hf1 = figure;
    set(hf1,'NumberTitle', 'off', ...
        'Name', sprintf('Participant: %s',ALLEEG(scounter).setname(1:s(1)-1)));
    scr_siz = get(0,'ScreenSize') ;
    pos = floor([scr_siz(3) scr_siz(4) scr_siz(3) scr_siz(4)]);
    set(hf1,'Position',pos,'Color',[1 1 1]);
    row = 8;
    clim = [-6 6];
    bline_t = T(blindx)
    
    for ecounter = 1:length(eindx)
    
        ax1(ecounter) = subplot(row,ceil(length(chans)/row),ecounter);
        imagesc('Parent',ax1(ecounter),'XData',T(tindx),'YData',foi,'CData',meanTFBL{1,ecounter},'CDataMapping','scaled')
        xlim([bline_t(1) Tnew(end)]);
        ax1(ecounter).YLim = [4 80];                    %Set Y-axis limites 
        ax1(ecounter).CLim = clim;                      %Set time-frequency power limits.
        ax1(ecounter).Layer = 'top';
        ax1(ecounter).XLabel.String = 'Time (seconds)'; %Set x-axis label
        ax1(ecounter).YLabel.String = 'Frequency (Hz)'; %Set y-axis label
        
        hold on
        plot(ax1(ecounter),T,coi,'w--','linewidth',2)
        title([chans{ecounter},' : ',bltype]);          
        colormap(ax1(ecounter),jet);
        cb(ecounter) = colorbar;                                 %Colorbar
        cb(ecounter).Label.String = ['ERSP (',bltype(1:2),' )']; %Colorbar title defined
        set(ax1(ecounter),'HitTest','on','SelectionHighlight','on','UserData',{T,tindx,foi,meanTFBL{1,ecounter},clim,chans{ecounter},bltype,coi},'Nextplot','replace');
        set(ax1(ecounter),'ButtonDownFcn',@plotsingle_tf)
    
    end

end

