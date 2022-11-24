
% Date: February 2014                       Programmer: Deirdre Bolger
% Function to carry out the basic EEG processing steps for a single EEGLAB
% dataset. It is possible to carry out the first 5 basic steps or to
% implement only an ICA on continuous data that has been cleaned.
% Inputs:
% - EEG : the EEGLAB data to be processed.
% - Params : a parameter structure with, at least, the following contents:
%                - srate : the sampling rate for downsampling
%                - fc_low : lower cutoff for filter
%                - fc_hi : high cutoff for filter
%                - electrodes : number of channels
% Output:
% -EEG : the output processed EEG dataset.
%*****************************************************************************************


prompt={'Preprocess/Chanrej/Epoch/ICA - choisir: ' 'Subject Numbers:'};
dlg_title='Choose Action';
num_lignes=[1;5];
action=inputdlg(prompt,dlg_title,num_lignes);


%% OPEN PARAMETERS FILE

fid=fopen(fullfile('Volumes','deepassport','Projects','Project-MotInter','Baseline=No_rs-filt-rref-ica','MotInter_EEG_LDT_param.txt'));      % il faut changer le chemin
mydata=textscan(fid,'%s %s');

for counter=1:length(mydata{1,1})                     % generate a parameters structure from the parameters text file
    Params.(genvarname(mydata{1,1}{counter}))=mydata{1,2}(counter);
end

if strcmp(action{1,1},'Preprocess')==1
    
    %% OPEN THE PARAMETERS TEXT FILE
    
    
    R=num2str(Params.references{1,1});
    if length(R)==2
        refs=str2double(R);
    elseif length(R)==4
        
        refs =[str2double(R(1:2)) str2double(R(3:4))];
    end
    
    dirdata=Params.Datadir{1,1};
    chaninfo='No';
    subs=cellstr(action{2,1});
    sujindx=cellfun(@str2double,subs);
    
    %% OPEN CURRENT EEGLAB SESSION
    
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;                %open eeglab session
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    
    %% START EEGLAB SESSION AND LOAD IN FILES - CAN TAKE ACCOUNT OF SUB-BLOCKS AND MERGE IF DESIRED
    
    for k=1:length(sujindx)
        
        %% VIEW ALL THE FILES INSIDE THIS FOLDER
        
        if sujindx(k)<10
            dir_open=strcat(dirdata,'Sujet',num2str(sujindx(k)),'\');                           %the directory of the current condition
        elseif sujindx(k)>=10
            dir_open=strcat(dirdata,'Sujet',num2str(sujindx(k)),'\');
        end
        display(dir_open);                                   %display the contents of the current folder
        allfiles= dir(dir_open);
        fileIndex = find(~[allfiles.isdir]);
        filenum=dir(strcat(dir_open,'*.bdf'));         %find all the *.bdf files in the current folder
        filenom={filenum.name};
        
        %% Determine the current directories
        
        if sujindx(k)<10
            currDir=dir_open;
            dir_save=strcat(dirdata,'Sujet',num2str(sujindx(k)),'\');
        elseif sujindx(k)>=10
            currDir=dir_open
            dir_save=strcat(dirdata,'Sujet',num2str(sujindx(k)),'\');
        end
        curr_title=filenom{1,k}(1:7);                                              %define the title of the current dataset
        
        if str2double(Params.blocknum{1,1})>1                             %if there are more than a single block i.e. there are sub-blocks
            
            
            for bloc_cnt=1:str2double(Params.blocknum{1,1})
                
                %bloc_ans=input('Do you wish to merge the blocks before preprocessing? (Yes/No)','s');
                if sujindx(k)<10
                    nf = strcat('Sujet', num2str(sujindx(k)),'liste',num2str(bloc_cnt),'.bdf');    % define the current file name.
                    titre = strcat('Sujet', num2str(sujindx(k)),'liste',num2str(bloc_cnt));       % Just the title of the file
                    titre_short=strcat('Sujet',num2str(sujindx(k)));
                elseif sujindx(k) >= 10
                    nf = strcat('Sujet', num2str(sujindx(k)),'liste',num2str(bloc_cnt),'.bdf');    % define the current file name.
                    titre = strcat('Sujet', num2str(sujindx(k)),'liste',num2str(bloc_cnt));
                    titre_short=strcat('Sujet',num2str(sujindx(k)));
                end
                
                fullDir = char(strcat(currDir, nf));     %Defines the full pathname
                
                if exist(fullDir,'file')==0     % if current file does not exist, continue to next iteration of i_bloc loop
                    
                    display('This file or repertoire does not exist!')
                    continue;
                end
                
                EEG = pop_biosig(fullDir, 'channels',[1:72], 'ref', [] ,'refoptions',{'keepref' 'off'} );
                [ALLEEG, EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(titre),'gui','off'); % Create a new dataset for the current raw datafile
                [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
                EEG = eeg_checkset( EEG );
                eeglab redraw
            end
            
            %% PREPARE INFORMATION TEXT-FILE
            fname=strcat(curr_title,'-info.txt');
            fdir=strcat(dir_save,fname);
            fid=fopen(fdir,'w');
            fprintf(fid,['---------',curr_title,'----------\n\n']);
            
            %% Merge the Lists Here
            
            EEG=pop_mergeset(ALLEEG,1:str2double(Params.blocknum{1,1}),0)
            EEG=eeg_checkset(EEG);
            MergeSet=strcat(titre_short,'_merged');
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(MergeSet),'gui','off'); % current set = xx;
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(MergeSet),'filepath',currDir);  % Saves a copy of the current resampled dataset to the current directory
            eeglab redraw
            
            
            %% Check the number of channels
            Enum=str2double(Params.electrodes{1,1});
            refs_string={EEG.chanlocs([refs]).labels};
            fprintf(fid,'Actual number of channels: %d\n',length(EEG.chanlocs));
            
            if length({EEG.chanlocs.labels})>Enum
                
                EEG = pop_select( EEG,'channel',1:Enum);
                [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(EEG.setname),'gui','off');
                EEG = eeg_checkset( EEG );
                EEG = pop_saveset( EEG, 'filename',char(EEG.setname),'filepath',currDir);  % Saves a copy of the current resampled dataset to the current directory
                eeglab redraw
                fprintf(fid,'Removed %d unnecessary channels\n',length(EEG.chanlocs)-Enum);
                
            end
            
            %% RESAMPLING OF DATA
            display('******************************Resampling from 20148Hz to 512Hz****************************');
            SR=str2double(Params.srate{1,1});
            fprintf(fid,'Downsampled from %fHz to %fHz\n',EEG.srate,SR);
            
            EEG = pop_resample(EEG, SR);   %resample the data at sampling rate defined, sr.
            EEG=eeg_checkset(EEG);
            Rs=strcat(MergeSet,'-rs');
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Rs),'gui','off'); % current set = xx;
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(Rs),'filepath',currDir);  % Saves a copy of the current resampled dataset to the current directory
            eeglab redraw
            
            
            %% FILTERING OF DATA USING WINDOWED SINC FILTER
            f_low=str2double(Params.fc_low{1,1});
            f_hi=str2double(Params.fc_hi{1,1});
            
            display('*********************************Bandpass filtering using a FIR windowed sinc filter***********************************')
            [M, dev]=pop_firwsord('blackman',SR, 2);
            [EEG,com,b]=pop_firws(EEG,'fcutoff',[f_low f_hi],'forder',M,'ftype','bandpass','wtype','blackman');
            fvtool(b);    % Visualise the filter characteristics
            Filtnom=strcat(Rs,'-filt');
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Filtnom),'gui','off');   %save the resampled data as a newdata set.
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(Filtnom),'filepath',currDir);
            eeglab redraw
            fprintf(fid,'Band-pass filtered %f - %fHz with %d-order fir sinc filter\n',f_low,f_hi,M);
            
            %% Reference the data
            display('***********************Rereference to Defined Channel:  does zero potential exist?*****************************')
            EEG=pop_reref(EEG, refs, 'method','standard','keepref','on');
            Rref=strcat(Filtnom,'-rref');
            EEG = eeg_checkset( EEG );
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Rref),'gui','off');   %save the resampled data as a newdata set.
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(Rref),'filepath',currDir);
            EEG = eeg_checkset( EEG );
            eeglab redraw
            fprintf(fid,['Rereferenced using channels ',EEG.chanlocs(refs(1)).labels,' & ',EEG.chanlocs(refs(2)).labels,'\n']);
            
            
            %% ADD CHANNEL INFORMATION
            chaninfo='Yes';
            if strcmp(chaninfo,'No')==1
                
                display('************Channel information already loaded!*************************')
                
            elseif strcmp(chaninfo,'Yes')==1
                display('*********************Adding Channel Information*****************************');
                chlocpath='C:\Users\bolger\Documents\MATLAB_tools\eeglab13_5_4b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp';
                EEG=pop_chanedit(EEG, 'lookup',chlocpath);            % Load channel path information
                [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                EEG = eeg_checkset( EEG );
                EEG = pop_saveset( EEG, 'filename',char(Rref),'filepath',currDir); %strcat('F:\BLRI\EEG\Projets_en_cours\Projet_MotInter\ExpEEG_Phase2\Data_Biosemi2\',groupe,'\All')
                eeglab redraw
                fprintf(fid,'Added channel information from %s\n',chlocpath);
            end
            
        elseif str2double(Params.blocknum{1,1})==1                                   %if there are no sub-blocks to load and merge
            
            
            fullDir = strcat(currDir,filenom{1,k});
            EEG = pop_biosig(fullDir, 'ref', [] ,'refoptions',{'keepref' 'off'});
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',curr_title,'gui','off'); % Create a new dataset for the current raw datafile
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = eeg_checkset( EEG );
            eeglab redraw
            
            %% PREPARE INFORMATION TEXT-FILE
            
            fname=strcat(curr_title,'-info.txt');
            fdir=strcat(dir_save,fname);
            fid=fopen(fdir,'w');
            fprintf(fid,['---------',curr_title,'----------\n\n']);
            
            %% CHECK THE NUMBER OF CHANNELS AND ADJUST IF NECESSARY
            
            Enum=str2double(Params.electrodes{1,1})
            fprintf(fid,'Actual number of channels: %d\n',length(EEG.chanlocs));
            
            if length({EEG.chanlocs.labels})>Enum
                EEG.setname=curr_title;
                EEG = pop_select( EEG,'channel',1:Enum);
                [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',curr_title,'gui','off');
                [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                EEG = eeg_checkset( EEG );
                EEG = pop_saveset( EEG, 'filename',curr_title,'filepath',dir_save);  % Saves a copy of the current resampled dataset to the current directory
                eeglab redraw
                
                fprintf(fid,'Removed %d unnecessary channels\n',length(EEG.chanlocs)-Enum);
                
            end
            
            %% RESAMPLING OF DATA
            display('******************************Resampling from 2048Hz to 512Hz****************************');
            SR=str2double(Params.srate{1,1});
            fprintf(fid,'Downsampled from %fHz to %fHz\n',EEG.srate,SR);
            
            EEG = pop_resample(EEG, SR);   %resample the data at sampling rate defined, sr.
            EEG=eeg_checkset(EEG);
            Rs=strcat(curr_title,'-rs');
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Rs),'gui','off'); % current set = xx;
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(Rs),'filepath',dir_save);  % Saves a copy of the current resampled dataset to the current directory
            eeglab redraw
            
            
            %% FILTERING OF DATA USING WINDOWED SINC FILTER
            f_low=str2double(Params.fc_low{1,1});
            f_hi=str2double(Params.fc_hi{1,1});
            
            display('*********************************Bandpass filtering using a FIR windowed sinc filter***********************************')
            [M, dev]=pop_firwsord('blackman',SR, 2);
            [EEG,com,b]=pop_firws(EEG,'fcutoff',[f_low f_hi],'forder',M,'ftype','bandpass','wtype','blackman');
            fvtool(b);    % Visualise the filter characteristics
            Filtnom=strcat(Rs,'-filt');
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Filtnom),'gui','off');   %save the resampled data as a newdata set.
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(Filtnom),'filepath',dir_save);
            eeglab redraw
            
            fprintf(fid,'Band-pass filtered %f - %fHz with %f-order fir win filter\n',f_low,f_hi,M);
            
            %% REREFERENCING DATA
            
            display('***********************Rereference to Average of Mastoids: does zero potential exist?*****************************')
            EEG=pop_reref(EEG, refs, 'method','standard');
            Rref=strcat(Filtnom,'_rref');
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Rref),'gui','off');   %save the resampled data as a newdata set.
            [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(Rref),'filepath',dir_save);
            eeglab redraw
            
            fprintf(fid,['Rereferenced using channels ',EEG.chanlocs(refs(1)).labels,' & ',EEG.chanlocs(refs(1)).labels,'\n']);
            
            %% ADD CHANNEL INFORMATION
            chaninfo='Yes';
            if strcmp(chaninfo,'No')==1
                
                display('************Channel information already loaded!*************************')
                
            elseif strcmp(chaninfo,'Yes')==1
                display('*********************Adding Channel Information*****************************');
                chlocpath='C:\Users\bolger\Documents\MATLAB_tools\eeglab13_5_4b\plugins\dipfit2.3\standard_BESA\standard-10-5-cap385.elp';
                EEG=pop_chanedit(EEG, 'lookup',chlocpath);            % Load channel path information
                [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                EEG = eeg_checkset( EEG );
                EEG = pop_saveset( EEG, 'filename',char(Rref),'filepath',dir_save);
                eeglab redraw
                fprintf(fid,'Added channel information from %s\n',chlocpath);
            end
            
        end
        
        %% Check that the triggers are in double format and not string format
        
        if ischar(EEG.event(5).type)
            disp('Trigger codes are in string format!')
            for trigi=1:length({EEG.event.type})
                EEG.event(trigi).type=str2double(EEG.event(trigi).type);
            end
        end
        beforeEp=CURRENTSET;
        %% SEGMENT THE DATA AND APPLY BASE-LINE CORRECTION
        %segment=input('Do you wish to epoch now? (yes/no)','s');
        segment='yes';
        
        if strcmp(segment,'yes')==1
            condnumber=str2double(Params.condnum{1,1});
            condcodes=zeros(condnumber,1);
            condtitle=cell(condnumber,1);
            allconds=unique([EEG.event.type]);                 %summary of all trigger codes found in the current dataset
            chans={EEG.chanlocs.labels};
            
            display('*********************Segmenting the Continuous Data*****************************');
            
            for cnt1=1:condnumber
                
                s1=strcat('condnom',num2str(cnt1));
                s2=strcat('condcode',num2str(cnt1));
                condtitle{cnt1,1}=Params.(genvarname(s1));
                condcodes(cnt1)=str2double(Params.(genvarname(s2)));
                H=find(ismember(allconds,condcodes(cnt1)));
                
                if isempty(H)==0
                    
                    display(strcat('*********Condition',num2str(condcodes(cnt1)),'****'));
                    
                    if cnt1>1
                        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',beforeEp,'study',0);
                        EEG = eeg_checkset( EEG );
                        eeglab redraw
                    end
                    
                    EEG.filepath=dir_save;
                    Enom=strcat(titre_short,'-',condtitle{cnt1,1})
                    EEG = pop_epoch( EEG, {condcodes(cnt1)}, [str2double(Params.wind_low{1,1}) str2double(Params.wind_hi{1,1})], 'newname', char(Enom), 'epochinfo', 'yes');
                    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Enom),'gui','off');
                    EEG = eeg_checkset( EEG );
                    EEG = pop_saveset( EEG, 'filename',char(Enom),'filepath',dir_save);
                    EEG = eeg_checkset( EEG );
                    
                    
                    disp('--------------------Baseline correction-----------------------------------');
                    Enom_bl=strcat(Enom,'-bl');
                    EEG = pop_rmbase( EEG, [Params.wind_low{1,1} 0]);
                    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Enom_bl),'gui','off');
                    EEG = eeg_checkset( EEG );
                    EEG = pop_saveset( EEG, 'filename',char(Enom_bl),'filepath',dir_save);
                    EEG = eeg_checkset( EEG );
                    eeglab redraw
                    
                    fprintf(fid,'\n*******************Condition %s:\n %d trials before cleaning\n',char(condtitle{cnt1,1}),size(EEG.data,3));
                    EEG = pop_autorej(EEG, 'nogui','on','threshold',75,'eegplot','on');
                    numrej=find(EEG.reject.rejauto);                 %number and indices of trials to reject
                    [EEG, eindx, measure,~] = pop_rejchan(EEG, 'elec',[1:64] ,'threshold',5,'norm','on','measure','kurt');
                    EEG.reject.rejkurtE=eindx;                          %indices of suggested electrodes to reject according to kurtosis.
                    
                    fprintf(fid,'Number of trials marked for rejection (>75mV): %d\n', length(numrej));
                    
                    
                    if ~isempty(eindx)
                        for cntr=1:length(eindx)
                            if cntr==1
                                fprintf(fid,'Bad electrodes according to kurtosis:  %s  ',chans{eindx(cntr)});
                            elseif cntr>1 && cntr<length(eindx)
                                fprintf(fid,' %s  ',chans{eindx(cntr)});
                            elseif cntr==length(eindx)
                                fprintf(fid,' %s \n ',chans{eindx(cntr)});
                            end
                            
                        end
                    end    %end if isempty if condition
                    
                elseif isempty(H)==1
                    
                    disp(strcat('***Condition ',num2str(condcodes(cnt1)),' not found in the current dataset!!'));
                    
                end
                
            end
            fprintf(fid,'\nLower and upper limits of baseline interval: %s s and %s s\n',char(Params.blc_low{1,1}),char(Params.blc_hi{1,1}));
            AllC=intersect(allconds,condcodes)';
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',beforeEp,'study',0);
            EEG = eeg_checkset( EEG );
            eeglab redraw
            Enom_all=strcat(titre_short,'-allconds');
            EEG = pop_epoch( EEG, {AllC(1,1) AllC(1,2) AllC(1,3) AllC(1,4) AllC(1,5) AllC(1,6) } ,[str2double(Params.wind_low{1,1}) str2double(Params.wind_hi{1,1})], 'newname', char(Enom_all), 'epochinfo', 'yes');
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Enom_all),'gui','off');
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(Enom_all),'filepath',dir_save);
            EEG = eeg_checkset( EEG );
            eeglab redraw
            
            %---------Baseline correction--------------------------------------
            Enom_all_bl=strcat(Enom_all,'-bl');
            EEG = pop_rmbase( EEG, [Params.wind_low{1,1} 0]);
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Enom_all_bl),'gui','off');
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(Enom_all_bl),'filepath',dir_save);
            EEG = eeg_checkset( EEG );
            eeglab redraw
        else
            disp('Leaving segmentation until later then------------------------');
        end
        
        fclose(fid);
        
    end
elseif strcmp(action{1,1},'Chanrej')==1
    %% REJECT CHANNELS
    currdir=EEG.filepath;
    prompt={'Channels to reject: ', 'Spherical Spline Interpolation (Y/N):'};
    dlg_title='Reject Channels';
    deflts = {'Fp1','N'};
    num_lignes=[10;3];
    chanrej_ans=inputdlg(prompt,dlg_title,num_lignes,deflts);
    options.resize='on';
    
    % Find the text file in the current dossier  (if it exists) and add the
    % titles of the channels to be rejected.
    dir_open=strcat(EEG.filepath,'\');                        %the directory of the current condition
    display(dir_open);                                   %display the contents of the current folder
    allfiles= dir(dir_open);
    fileIndex = find(~[allfiles.isdir]);
    filenum=dir(strcat(dir_open,'*.txt'));
    fid=fopen(strcat(dir_open,filenum.name),'a');
    fprintf(fid,' **********************Channels Rejected*********************************\n');
    chanindx = zeros(size(chanrej_ans{1,1},1),1);   %for the channel indices
    for c=1:size(chanrej_ans{1,1},1)
        
        if c==size(chanrej_ans{1,1},1)
            fprintf(fid,'%s\n\n',chanrej_ans{1,1}(c,:));
            a=chanrej_ans{1,1}(c,:);
            chanindx(c)=find(strcmp({EEG.chanlocs.labels},cellstr(a)));
            clear a;
        elseif c<size(chanrej_ans{1,1},1)
            fprintf(fid,'%s\t',chanrej_ans{1,1}(c,:));
            a=chanrej_ans{1,1}(c,:);
            chanindx(c)=find(strcmp({EEG.chanlocs.labels},cellstr(a)));  %{EEG.chanlocs.labels}
        end
    end
    
    here=CURRENTSET;
    title_rej = strcat(EEG.setname,'-chanrej');
    EEG=pop_select(EEG,'nochannel',chanindx);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(title_rej),'gui','off'); % current set = xx;
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset( EEG, 'filename',char(title_rej),'filepath',currdir);  % Saves a copy of the current resampled dataset to the current directory
    eeglab redraw
    
    if strcmp(chanrej_ans{2,1},'Y')==1
        disp('***************************Carrying out Spherical Spline Interpolation on current dataset***************************');
        
        spline_note='Carried out spherical spline interpolation';
        fprintf(fid,'%s\n\n',spline_note);
        fclose(fid);
        
        title_inter=strcat(EEG.setname,'-ssinterp');
        EEG = pop_interp(EEG, ALLEEG(here).chanlocs, 'spherical');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(title_inter),'gui','off'); % current set = xx;
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename',char(title_inter),'filepath',currdir);  % Saves a copy of the current resampled dataset to the current directory
        eeglab redraw
    else
        disp('***************No Interpolation for the moment***************************');
    end
    
elseif strcmp(action{1,1},'Epoch')==1
    
    %% Check that the triggers are in double format and not string format
    
    if ischar(EEG.event(5).type)
        disp('Trigger codes are in string format!')
        for trigi=1:length({EEG.event.type})
            EEG.event(trigi).type=str2double(EEG.event(trigi).type);
        end
    end
    
    bl_ans=input('Baseline correction?(Y/N)','s');
    nombase=EEG.setname(1:9);
    here = CURRENTSET;
    curdir=EEG.filepath;
    %curdir=strcat(curdir,filesep);
    condcodes=zeros(str2double(Params.condnum),1);
    condnames=cell(str2double(Params.condnum),1);
    for cindx=1:str2double(Params.condnum)
        condcodes(cindx)=str2double(Params.(genvarname(strcat('condcode',num2str(cindx)))));
        condnames{cindx,1}=Params.(genvarname(strcat('condnom',num2str(cindx))));
    end
    condcurr=unique([EEG.event.type]);
    condcurr=condcurr([~isnan(condcurr)]);
    condcurr = intersect(condcodes,condcurr);
    
    for cindx2=1:length(condcurr)
        
        if cindx2>1
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',here,'study',0);
            EEG = eeg_checkset( EEG );
            eeglab redraw
        end
        
           nomnew=strcat(nombase,'-',condnames{ismember(condcodes, condcurr(cindx2)),1});
        EEG = pop_epoch( EEG, {condcurr(cindx2)}, [str2double(Params.wind_low{1,1}) str2double(Params.wind_hi{1,1})], 'newname',char(nomnew) ,'epochinfo', 'yes');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(nomnew),'gui','off');
        EEG.condition=char(condnames{ismember(condcodes, condcurr(cindx2)),1});
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename',char(nomnew),'filepath',curdir);
        EEG = eeg_checkset( EEG );
        eeglab redraw
        
        %---------Baseline correction--------------------------------------
        if strcmp(bl_ans,'Y')
            nomnew_bl=strcat(nomnew,'-bl');
            EEG = pop_rmbase( EEG, [Params.blc_low{1,1} 0]);
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(nomnew_bl),'gui','off');
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(nomnew_bl),'filepath',curdir);
            EEG = eeg_checkset( EEG );
            eeglab redraw
        else
            disp('-----------Skipping baseline correction--------------------');
        end
        
    end
    
    seg_ans=input('Segment all conditions together?(Y/N)','s');
    if strcmp(seg_ans,'Y')
        disp('----------Segmentation all conditions together--------------')
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',here,'study',0);
        EEG = eeg_checkset( EEG );
        eeglab redraw
        
        Enom_all=strcat(nombase,'-allconds');
        EEG = pop_epoch( EEG, num2cell(condcurr) ,[str2double(Params.wind_low{1,1}) str2double(Params.wind_hi{1,1})], 'newname', char(Enom_all), 'epochinfo', 'yes');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Enom_all),'gui','off');
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename',char(Enom_all),'filepath',curdir);
        EEG = eeg_checkset( EEG );
        eeglab redraw
        
        %---------Baseline correction--------------------------------------
        if strcmp(bl_ans,'Y')
            Enom_all_bl=strcat(Enom_all,'-bl');
            EEG = pop_rmbase( EEG, [Params.blc_low{1,1} 0]);
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',char(Enom_all_bl),'gui','off');
            EEG = eeg_checkset( EEG );
            EEG = pop_saveset( EEG, 'filename',char(Enom_all_bl),'filepath',curdir);
            EEG = eeg_checkset( EEG );
            eeglab redraw
        else
            disp('-----------Skipping baseline correction--------------------');
        end
    end
    
    
    
elseif strcmp(action{1,1},'ICA')==1
    
    %% CARRY OUT ICA
    
    
    for datacnt=1:length(ALLEEG)
        EEG=ALLEEG(datacnt);
        icalen=length(EEG.chanlocs)-8;
        display('*********************Carrying out ica: patience!****************************');
        [weights, sphere,compvars,bias,signs,lrates, activations]=runica(EEG.data(1:icalen,:), 'extended',0);
        
        EEG.icaweights=weights;
        EEG.icasphere=sphere;
        EEG.icachansind=1:icalen;
        EEG.icaact=activations;
        [icaprojdata]=icaproj(EEG.data(1:icalen,:),weights,1:icalen);
        EEG.icawinv=inv(weights*sphere);
        
        [ALLEEG EEG]=eeg_store(ALLEEG, EEG, CURRENTSET);
        EEG=eeg_checkset(EEG)
        reply='yes';
        
        if strcmp(reply,'yes')==1
            EEG=pop_saveset(EEG, 'filename',char(strcat(char(EEG.setname),'-ica')),'filepath',EEG.filepath); %change EEG.filepath
        end
        eeglab redraw
        EEG=[];
    end
    
end
