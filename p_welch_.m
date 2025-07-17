%%  Welch power spectral density estimate

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  12/2023
% Last update: 01/2024


% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m


%%
fprintf('\n Welch power spectral density estimate ... \n');

%% pwelch
%  [pxx,f] = pwelch(x,window,noverlap,f,fs)

% - cell
% - 1st column    - > baseline
% - 2nd column    - > cs events

% in each cell
% - rows          - > Hz
% - columns       - > channels
% 3th dimentions  - > CS-Trials

pw = [];

% Time window
pw.full_trial.parameters.baseline_timewin    = 2000; % in ms
pw.full_trial.parameters.STIM_Trials_timewin = 2000; % in ms
pw.full_trial.parameters.ITI_timewin         = 2000; % in ms

% Convert time window to points
pw.full_trial.parameters.baseline_timewinpnts     = hamming(round(pw.full_trial.parameters.baseline_timewin/(1000/parameters.decimated_srate)));
pw.full_trial.parameters.STIM_Trials_timewinpnts  = hamming(round(pw.full_trial.parameters.STIM_Trials_timewin/(1000/parameters.decimated_srate)));
pw.full_trial.parameters.ITI_timewinpnts          = hamming(round(pw.full_trial.parameters.ITI_timewin/(1000/parameters.decimated_srate)));

% nFFT
pw.full_trial.parameters.nFFT = 2^15; %4096; %2^nextpow2(pw.full_trial.parameters.baseline_timewinpnts));

% Number of overlap samples
pw.full_trial.parameters.overlap = 90;
pw.full_trial.parameters.baseline_noverlap = floor(pw.full_trial.parameters.overlap*0.01 * pw.full_trial.parameters.baseline_timewin);
pw.full_trial.parameters.STIM_Trials_noverlap = floor(pw.full_trial.parameters.overlap*0.01 * pw.full_trial.parameters.STIM_Trials_timewin);
pw.full_trial.parameters.ITI_noverlap = floor(pw.full_trial.parameters.overlap*0.01 * pw.full_trial.parameters.ITI_timewin);


% Baseline
not1 = 6;

for ii = 1:size(data.lfp{not1, 1},1)

    if ii == 1
        [pw.full_trial.Pxx{1,1}(ii,:),pw.full_trial.freq_baseline] = pwelch(data.lfp{not1, 1}(ii,B_clean{ms}(1):B_clean{ms}(2)),pw.full_trial.parameters.baseline_timewinpnts,pw.full_trial.parameters.baseline_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

    else
        pw.full_trial.Pxx{1,1}(ii,:) = pwelch(data.lfp{not1, 1}(ii,B_clean{ms}(1):B_clean{ms}(2)),pw.full_trial.parameters.baseline_timewinpnts,pw.full_trial.parameters.baseline_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

    end

end

clear ('ii')


% 4Hz or 8Hz opto stimulation
% Choose data from data.lfp according pre_processing.m define
not2 = 7;

for jj = 1:size(data.lfp{not2, 1},3)
    for ii = 1:size(data.lfp{not2, 1},1)
        if ii == 1
            [pw.full_trial.Pxx{2,1}(ii,:,jj),pw.full_trial.freq_STIM_Trials] = pwelch(data.lfp{not2, 1}(ii,:,jj),pw.full_trial.parameters.STIM_Trials_timewinpnts,pw.full_trial.parameters.STIM_Trials_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{2,1}(ii,:,jj) = pwelch(data.lfp{not2, 1}(ii,:,jj),pw.full_trial.parameters.STIM_Trials_timewinpnts,pw.full_trial.parameters.STIM_Trials_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        end

    end
end

clear ('jj','ii')

% % ITI or non-freezing epochs
% Choose data from data.lfp according pre_processing.m define
not3 = 8;

for jj = 1:size(data.lfp{not3, 1},3)
    for ii = 1:size(data.lfp{not3, 1},1)
        if ii == 1
            [pw.full_trial.Pxx{3,1}(ii,:,jj),pw.full_trial.freq_ITI] = pwelch(data.lfp{not3, 1}(ii,:,jj),pw.full_trial.parameters.ITI_timewinpnts,pw.full_trial.parameters.ITI_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        else
            pw.full_trial.Pxx{3,1}(ii,:,jj) = pwelch(data.lfp{not3, 1}(ii,:,jj),pw.full_trial.parameters.ITI_timewinpnts,pw.full_trial.parameters.ITI_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);

        end

    end
end

% % 4Hz or 8Hz opto stimulation
% % Choose data from data.lfp according pre_processing.m define
% not4 = 9;
% 
% for jj = 1:size(data.lfp{not4, 1},3)
%     for ii = 1:size(data.lfp{not4, 1},1)
%         if ii == 1
%             [pw.full_trial.Pxx{4,1}(ii,:,jj),pw.full_trial.freq_STIM_Trials] = pwelch(data.lfp{not4, 1}(ii,:,jj),pw.full_trial.parameters.STIM_Trials_timewinpnts,pw.full_trial.parameters.STIM_Trials_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);
% 
%         else
%             pw.full_trial.Pxx{4,1}(ii,:,jj) = pwelch(data.lfp{not4, 1}(ii,:,jj),pw.full_trial.parameters.STIM_Trials_timewinpnts,pw.full_trial.parameters.STIM_Trials_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);
% 
%         end
% 
%     end
% end
% 
% clear ('jj','ii')
% 
% % % ITI or non-freezing epochs
% % Choose data from data.lfp according pre_processing.m define
% not5 = 10;
% 
% for jj = 1:size(data.lfp{not5, 1},3)
%     for ii = 1:size(data.lfp{not5, 1},1)
%         if ii == 1
%             [pw.full_trial.Pxx{5,1}(ii,:,jj),pw.full_trial.freq_ITI] = pwelch(data.lfp{not5, 1}(ii,:,jj),pw.full_trial.parameters.ITI_timewinpnts,pw.full_trial.parameters.ITI_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);
% 
%         else
%             pw.full_trial.Pxx{5,1}(ii,:,jj) = pwelch(data.lfp{not5, 1}(ii,:,jj),pw.full_trial.parameters.ITI_timewinpnts,pw.full_trial.parameters.ITI_noverlap,pw.full_trial.parameters.nFFT,parameters.decimated_srate);
% 
%         end
% 
%     end
% end

clear ('jj','ii')

% figure
% plot(pw.full_trial.freq_baseline,pw.full_trial.Pxx{1,1}(1,:))
% hold on
% plot(pw.full_trial.freq_STIM_Trials,mean(pw.full_trial.Pxx{2,1}(1,:,1:3),3))
% plot(pw.full_trial.freq_STIM_Trials,mean(pw.full_trial.Pxx{3,1}(1,:,1:3),3))
%
% xlim([3 12])

clear ('not1','not2','not3')

%% Choosing Channels and averaging

PL_  = [1];  % mPFC PL
IL_  = [3];  % mPFC IL
dHPC = [15]; % dHPC


for ii = 1:size(pw.full_trial.Pxx,1)

    pw.full_trial.Pxx{ii,2}(1,:,:) = mean(pw.full_trial.Pxx{ii,1}(PL_,:,:),1);
    pw.full_trial.Pxx{ii,2}(2,:,:) = mean(pw.full_trial.Pxx{ii,1}(IL_,:,:),1);
    pw.full_trial.Pxx{ii,2}(3,:,:) = mean(pw.full_trial.Pxx{ii,1}(dHPC,:,:),1);

    if ii == 1
        pw.full_trial.Pxx{ii,2} = squeeze(pw.full_trial.Pxx{ii,2});
    end
end

clear('PL_', 'IL_', 'dHPC', 'ii')

%% Normalization and trials average -> 2 - 12 HZ

pw.full_trial.Pxx_TotalPower_norm_2_12_Hz = [];
pw.full_trial.Pxx_baseline_norm_2_12_Hz   = [];
pw.full_trial.Pxx_z_norm_2_12_Hz          = [];

% Frequency resolution. According to the fft time window
steps                                     = diff(pw.full_trial.freq_baseline);

% Set freuency range
pw.full_trial.parameters.frex_2_12_Hz      = 2:steps(1):12;
pw.full_trial.parameters.frex_idx_2_12_Hz  = dsearchn(pw.full_trial.freq_baseline,pw.full_trial.parameters.frex_2_12_Hz');
pw.full_trial.parameters_freq_v_2_12_      = pw.full_trial.freq_baseline(pw.full_trial.parameters.frex_idx_2_12_Hz);


% Plexon rescale to uV
% The original files *.nex5 did not have the correct scale.
scale = 1/0.153;

% Total power
for jj = 1:size(pw.full_trial.Pxx,1)
    for ii = 1:size(pw.full_trial.Pxx,2)

        % Total power normalization
        pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{jj,ii} = scale.*(pw.full_trial.Pxx{jj,ii}(:,pw.full_trial.parameters.frex_idx_2_12_Hz,:)./sum(pw.full_trial.Pxx{jj,ii}(:,pw.full_trial.parameters.frex_idx_2_12_Hz,:),2));
        
        % Baseline Normalization
        pw.full_trial.Pxx_baseline_norm_2_12_Hz{jj,ii} = pw.full_trial.Pxx{jj,ii}(:,pw.full_trial.parameters.frex_idx_2_12_Hz,:) ./ pw.full_trial.Pxx{1,ii}(:,pw.full_trial.parameters.frex_idx_2_12_Hz,:);
        %pw.full_trial.Pxx_baseline_norm_2_12_Hz{jj,ii} = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{jj,ii}./ pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,ii};
        
        % Zscore Normalization
        % pw.full_trial.Pxx_z_norm_2_12_Hz{jj,ii} = ((pw.full_trial.Pxx{jj,ii}(:,pw.full_trial.parameters.frex_idx_2_12Hz,:) - mean(pw.full_trial.Pxx{1,ii}(:,pw.full_trial.parameters.frex_idx_2_12Hz,:),2))./...
        %     std(pw.full_trial.Pxx{1,ii}(:,pw.full_trial.parameters.frex_idx_2_12Hz,:),[],2)); %zscore

        pw.full_trial.Pxx_z_norm_2_12_Hz{jj,ii} = ((pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{jj,ii} - mean(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,ii},2))./...
            std(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,ii},[],2));
        
    end
end



%%
% on/off normalization -> 
% The power spectrum during laser on and off were normalized to the maxima and minima of laser off data
% Jordan S. Farrell et al. 
% Science 374, 1492 (2021)
% DOI: 10.1126/science.abh4272


pw.full_trial.Pxx_on_off_2_12_Hz = [];

        pw.full_trial.Pxx_on_off_2_12_Hz{1,1} = (pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{2,1} - min(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,1},[],2)) ./ (max(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,1},[],2) - min(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,1},[],2));
        pw.full_trial.Pxx_on_off_2_12_Hz{1,2} = (pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{2,2} - min(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,2},[],2)) ./ (max(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,2},[],2) - min(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,2},[],2));

        pw.full_trial.Pxx_on_off_2_12_Hz{2,1} = (pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{3,1} - min(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,1},[],2)) ./ (max(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,1},[],2) - min(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,1},[],2));
        pw.full_trial.Pxx_on_off_2_12_Hz{2,2} = (pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{3,2} - min(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,2},[],2)) ./ (max(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,2},[],2) - min(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,2},[],2));






clear('ii','jj','scale','steps')

%% Plot to check 

data_2_plot      = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{2, 1};
%data_2_plot       = pw.full_trial.Pxx_baseline_norm_2_12_Hz{2, 1};
%data_2_plot      = pw.full_trial.Pxx_z_norm_2_12_Hz{2, 1};


figure;
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','w');
axis('square')
sgtitle('\fontsize{18} \bf Extinction')



for ii = 1:size(data_2_plot,1)

    subplot(4,8,ii);
    %subplot(1,3,ii);
    hold on
    plot(pw.full_trial.parameters.frex_2_12_Hz,mean(data_2_plot(ii,:,5),3),'linewidth',4,'color',[.6 .6 .6] )

    %xlim([2 12])
    %ylim([-3 20])
    %ylim([0 .1])
    %ylim([0 3])
    ylabel({'Normalized Amplitude (a.u.)';[]})
    xlabel('Hz')
    %yticks([0 0.05 0.1 0.15 0.2])
    box off
    title({'\fontsize{16} \bf Channel: ' num2str(ii)})

end

clear ('ii',"data_2_plot")


%% Check Trials

data_2_plot_1      = pw.full_trial.Pxx_on_off_2_12_Hz;

figure;
%set(gcf, 'Position', get(0, 'Screensize'));
set(gcf, 'Position', [1683, -377, 1249,1068]);
set(gcf,'color','w');
axis('square')
%sgtitle({'\fontsize{18} \bf Habituation',[]})
sgtitle('\fontsize{18} \bf Extinction')



for ii = 1:45

    subplot(5,9,ii);
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{1,2}(2,:,ii),'linewidth',3,'color',[.2, .6, 1])
    hold on
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,2}(2,:,ii),'linewidth',3,'color',[.7 .7 .7])
    %set(gca, 'YTick', [],'YColor','none')
    box off
    xlim([3 12])
    ylim([0 1])
    title({['\fontsize{14}Trial ' num2str(ii)],[]})
end


clear('data_2_plot_1','ii')

%% Remove noise trials

nan_cell_1 = nan(size(pw.full_trial.Pxx_z_norm_2_12_Hz{1, 1}));
nan_cell_2 = nan(size(pw.full_trial.Pxx_z_norm_2_12_Hz{1, 2}));


noise_trial_ = [2 4 5 6 11 13];

for ii = 2:size(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz,1)
    for jj = 1:size(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz,2)
        for tt = 1:size(noise_trial_,2)

            if jj == 1
                pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{ii,jj}(:,:,noise_trial_(1,tt)) = nan_cell_1;
                pw.full_trial.Pxx_baseline_norm_2_12_Hz{ii,jj}(:,:,noise_trial_(1,tt))   = nan_cell_1;
                pw.full_trial.Pxx_z_norm_2_12_Hz{ii,jj}(:,:,noise_trial_(1,tt))          = nan_cell_1;
            else
                pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{ii,jj}(:,:,noise_trial_(1,tt)) = nan_cell_2;
                pw.full_trial.Pxx_baseline_norm_2_12_Hz{ii,jj}(:,:,noise_trial_(1,tt))   = nan_cell_2;
                pw.full_trial.Pxx_z_norm_2_12_Hz{ii,jj}(:,:,noise_trial_(1,tt))          = nan_cell_2;
            end

        end
    end
end

clear('nan_cell_1','nan_cell_2','noise_trial_','ii','jj','tt')

%%

nan_cell_2 = nan(1,328);
noise_trial_ = [2 4 5 6 7 9 10 11 13 15 16 17 18 19 21 22 24 25 27 28 29 30 31 32 36 37 38 42 44 45];

for ii = 2:size(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz,1)
    for jj = 2
        for tt = 1:size(noise_trial_,2)
            pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{ii,jj}(3,:,noise_trial_(1,tt)) = nan_cell_2;
            pw.full_trial.Pxx_baseline_norm_2_12_Hz{ii,jj}(3,:,noise_trial_(1,tt))   = nan_cell_2;
            pw.full_trial.Pxx_z_norm_2_12_Hz{ii,jj}(3,:,noise_trial_(1,tt))          = nan_cell_2;
        end
    end
end


%% Averaging 5 trials Block

blk_ = 1:5:size(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{2,2},3)+1;


for ii = 2:size(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz,1)
        for tt = 1:size(blk_,2)-1
          
            pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{ii,3}(:,:,tt) = mean(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{ii,2}(:,:,blk_(tt):blk_(tt+1)-1),3,"omitnan");
            pw.full_trial.Pxx_baseline_norm_2_12_Hz{ii,3}(:,:,tt)   = mean(pw.full_trial.Pxx_baseline_norm_2_12_Hz{ii,2}(:,:,blk_(tt):blk_(tt+1)-1),3,"omitnan");
            pw.full_trial.Pxx_z_norm_2_12_Hz{ii,3}(:,:,tt)          = mean(pw.full_trial.Pxx_z_norm_2_12_Hz{ii,2}(:,:,blk_(tt):blk_(tt+1)-1),3,"omitnan");

        end
end


for tt = 1:size(blk_,2)-1
    pw.full_trial.Pxx_on_off_2_12_Hz{1,3}(:,:,tt) = mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,2}(:,:,blk_(tt):blk_(tt+1)-1),3,"omitnan");
    pw.full_trial.Pxx_on_off_2_12_Hz{2,3}(:,:,tt) = mean(pw.full_trial.Pxx_on_off_2_12_Hz{2,2}(:,:,blk_(tt):blk_(tt+1)-1),3,"omitnan");
end

%% Plot to check - 5 trials Block - on/off

data_2_plot_1      = pw.full_trial.Pxx_on_off_2_12_Hz;

figure;
%set(gcf, 'Position', get(0, 'Screensize'));
set(gcf, 'Position', [1683, -377, 1822,1068]);
set(gcf,'color','w');
axis('square')
%sgtitle({'\fontsize{18} \bf Habituation',[]})
sgtitle('\fontsize{18} \bf Extinction')



% mPFC IL
for ii = 1:9

    subplot(3,10,ii);
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{1,3}(1,:,ii),'linewidth',3,'color',[.2, .6, 1])
    hold on
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,3}(1,:,ii),'linewidth',3,'color',[.7 .7 .7])
    %set(gca, 'YTick', [],'YColor','none')
    box off
    xlim([3 12])
    %ylim([0 .15])
    title({['\fontsize{14}CS-Block ' num2str(ii)],[]})
end

subplot(3,10,10);
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{1,3}(1,:,:),3,"omitnan"),'linewidth',3,'color',[.2, .6, 1])
hold on
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,3}(1,:,:),3,"omitnan"),'linewidth',3,'color',[.7 .7 .7])
%set(gca, 'YTick', [])
xlim([3 12])
ylim([-0.05 3])
title({'\fontsize{14}Averaged Trials',[]})




% mPFC PL
for ii = 1:9

    subplot(3,10,ii+10);
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{1,3}(2,:,ii),'linewidth',3,'color',[.2, .6, 1])
    hold on
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,3}(2,:,ii),'linewidth',3,'color',[.7 .7 .7])
    %set(gca, 'YTick', [],'YColor','none')
    box off
    xlim([3 12])
    %ylim([0 .15])
end

subplot(3,10,20);
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{1,3}(2,:,:),3,"omitnan"),'linewidth',3,'color',[.2, .6, 1])
hold on
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,3}(2,:,:),3,"omitnan"),'linewidth',3,'color',[.7 .7 .7])
%set(gca, 'YTick', [])
xlim([3 12])
ylim([-0.05 3])




% dHPC

for ii = 1:9

    subplot(3,10,ii+20);
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{1,3}(3,:,ii),'linewidth',3,'color',[.2, .6, 1])
    hold on
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,3}(3,:,ii),'linewidth',3,'color',[.7 .7 .7])
    %set(gca, 'YTick', [],'YColor','none')
    box off
    xlim([3 12])
    %ylim([0 .15])
end

subplot(3,10,30);
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{1,3}(3,:,:),3,"omitnan"),'linewidth',3,'color',[.2, .6, 1])
hold on
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,3}(3,:,:),3,"omitnan"),'linewidth',3,'color',[.7 .7 .7])
%set(gca, 'YTick', [])
xlim([3 12])
ylim([-0.05 3])
legend('CS + 8 HZ Stim','ITI')
legend('boxoff')



clear ('ii',"data_2_plot_2","data_2_plot_3","data_2_plot_1")

%% Plot to check - 5 trials Block 

data_2_plot_1      = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz;
% data_2_plot_1      = pw.full_trial.Pxx_baseline_norm_2_12_Hz;
% data_2_plot_1      = pw.full_trial.Pxx_z_norm_2_12_Hz;


figure;
%set(gcf, 'Position', get(0, 'Screensize'));
set(gcf, 'Position', [1683, -377, 1822,1068]);
set(gcf,'color','w');
axis('square')
%sgtitle({'\fontsize{18} \bf Habituation',[]})
sgtitle('\fontsize{18} \bf Extinction')



% mPFC IL
subplot(3,11,1);
plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{1,2}(1,:),'linewidth',2,'color','k')
ylabel({'\fontsize{16}mPFC IL';'Normalized Amplitude (a.u.)';[]})
box off
xlim([3 12])
ylim([0 .1])
title({'\fontsize{14}Baseline',[]})

for ii = 1:9

    subplot(3,11,ii+1);
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,3}(1,:,ii),'linewidth',3,'color',[.2, .6, 1])
    hold on
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{3,3}(1,:,ii),'linewidth',3,'color',[.7 .7 .7])
    set(gca, 'YTick', [],'YColor','none')
    box off
    xlim([3 12])
    ylim([0 .1])
    title({['\fontsize{14}CS-Block ' num2str(ii)],[]})
end

subplot(3,11,11);
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,3}(1,:,:),3),'linewidth',3,'color',[.2, .6, 1])
hold on
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{3,3}(1,:,:),3),'linewidth',3,'color',[.7 .7 .7])
set(gca, 'YTick', [])
xlim([3 12])
ylim([0 .1])
title({'\fontsize{14}Averaged Trials',[]})




% mPFC PL
subplot(3,11,12);
plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{1,2}(2,:),'linewidth',2,'color','k')
ylabel({'\fontsize{16}mPFC PL';'Normalized Amplitude (a.u.)';[]})
box off
xlim([3 12])
ylim([0 .1])

for ii = 1:9

    subplot(3,11,ii+12);
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,3}(2,:,ii),'linewidth',3,'color',[.2, .6, 1])
    hold on
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{3,3}(2,:,ii),'linewidth',3,'color',[.7 .7 .7])
    set(gca, 'YTick', [],'YColor','none')
    box off
    xlim([3 12])
    ylim([0 .1])
end

subplot(3,11,22);
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,3}(2,:,:),3),'linewidth',3,'color',[.2, .6, 1])
hold on
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{3,3}(2,:,:),3),'linewidth',3,'color',[.7 .7 .7])
set(gca, 'YTick', [])
xlim([3 12])
ylim([0 .1])




% dHPC
subplot(3,11,23);
plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{1,2}(3,:),'linewidth',2,'color','k')
ylabel({'\fontsize{16}dHPC';'Normalized Amplitude (a.u.)';[]})
box off
xlim([3 12])
ylim([0 .15])

for ii = 1:9

    subplot(3,11,ii+23);
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,3}(3,:,ii),'linewidth',3,'color',[.2, .6, 1])
    hold on
    plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{3,3}(3,:,ii),'linewidth',3,'color',[.7 .7 .7])
    set(gca, 'YTick', [],'YColor','none')
    box off
    xlim([3 12])
    ylim([0 .15])
end

subplot(3,11,33);
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,3}(3,:,:),3),'linewidth',3,'color',[.2, .6, 1])
hold on
plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{3,3}(3,:,:),3),'linewidth',3,'color',[.7 .7 .7])
set(gca, 'YTick', [])
xlim([3 12])
ylim([0 .15])
legend('CS + 8 HZ Stim','ITI')
legend('boxoff')



clear ('ii',"data_2_plot_2","data_2_plot_3","data_2_plot_1")

%% Plot to check - 10 first and 10 last trials

% % data_2_plot_1      = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz;
% % data_2_plot_2      = pw.full_trial.Pxx_baseline_norm_2_12_Hz;
% % data_2_plot_3      = pw.full_trial.Pxx_z_norm_2_12_Hz;
% 
% best_trials1 = 21:30;
% best_trials2 = 36:45;
% 
% figure;
% %set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf, 'Position', [1683, -377, 1249,1068]);
% set(gcf,'color','w');
% axis('square')
% sgtitle({'\fontsize{18} \bf Habituation',[]})
% %sgtitle('\fontsize{18} \bf Extinction')
% 
% 
% 
% % mPFC IL
% subplot(6,12,[1 13]);
% plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{1,2}(1,:),'linewidth',2,'color','k')
% ylabel({'\fontsize{16}mPFC IL';'Normalized Amplitude (a.u.)';[]})
% box off
% xlim([3 12])
% ylim([0 .15])
% title({'\fontsize{14}Baseline',[]})
% 
% for ii = 1:10
% 
%     subplot(6,12,ii+1);
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,2}(1,:,ii),'linewidth',3,'color',[.2, .6, 1])
%     hold on
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{3,2}(1,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     ylim([0 .15])
%     title({['\fontsize{14}Trial ' num2str(ii)],[]})
% end
% 
% for ii = 36:45
% 
%     subplot(6,12,ii-22);
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,2}(1,:,ii),'linewidth',3,'color',[.2, .6, 1])
%     hold on
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{3,2}(1,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     ylim([0 .15])
%     title({['\fontsize{14}Trial ' num2str(ii)],[]})
% end
% 
% subplot(6,12,12);
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,2}(1,:,best_trials1),3),'linewidth',3,'color',[.2, .6, 1])
% hold on
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{3,2}(1,:,best_trials1),3),'linewidth',3,'color',[.7 .7 .7])
% set(gca, 'YTick', [])
% xlim([3 12])
% ylim([0 .15])
% title({'\fontsize{14}Averaged Best Trials',[]})
% 
% subplot(6,12,24);
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,2}(1,:,best_trials2),3),'linewidth',3,'color',[.2, .6, 1])
% hold on
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{3,2}(1,:,best_trials2),3),'linewidth',3,'color',[.7 .7 .7])
% set(gca, 'YTick', [])
% xlim([3 12])
% ylim([0 .15])
% 
% 
% 
% % mPFC PL
% subplot(6,12,[25 37]);
% plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{1,2}(2,:),'linewidth',2,'color','k')
% ylabel({'\fontsize{16}mPFC PL';'Normalized Amplitude (a.u.)';[]})
% box off
% xlim([3 12])
% ylim([0 .15])
% 
% for ii = 1:10
% 
%     subplot(6,12,ii+25);
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,2}(2,:,ii),'linewidth',3,'color',[.2, .6, 1])
%     hold on
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{3,2}(2,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     ylim([0 .15])
%     title({['\fontsize{14}Trial ' num2str(ii)],[]})
% end
% 
% for ii = 36:45
% 
%     subplot(6,12,ii+2);
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,2}(2,:,ii),'linewidth',3,'color',[.2, .6, 1])
%     hold on
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{3,2}(2,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     ylim([0 .15])
%     title({['\fontsize{14}Trial ' num2str(ii)],[]})
% end
% 
% subplot(6,12,36);
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,2}(2,:,best_trials1),3),'linewidth',3,'color',[.2, .6, 1])
% hold on
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{3,2}(2,:,best_trials1),3),'linewidth',3,'color',[.7 .7 .7])
% set(gca, 'YTick', [])
% xlim([3 12])
% ylim([0 .15])
% 
% subplot(6,12,48);
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,2}(2,:,best_trials2),3),'linewidth',3,'color',[.2, .6, 1])
% hold on
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{3,2}(2,:,best_trials2),3),'linewidth',3,'color',[.7 .7 .7])
% set(gca, 'YTick', [])
% xlim([3 12])
% ylim([0 .15])
% 
% 
% 
% % dHPC
% subplot(6,12,[49 61]);
% plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{1,2}(3,:),'linewidth',2,'color','k')
% ylabel({'\fontsize{16}dHPC';'Normalized Amplitude (a.u.)';[]})
% box off
% xlim([3 12])
% ylim([0 .15])
% 
% for ii = 1:10
% 
%     subplot(6,12,ii+49);
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,2}(3,:,ii),'linewidth',3,'color',[.2, .6, 1])
%     hold on
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{3,2}(3,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     ylim([0 .15])
%     title({['\fontsize{14}Trial ' num2str(ii)],[]})
% end
% 
% for ii = 36:45
% 
%     subplot(6,12,ii+26);
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{2,2}(3,:,ii),'linewidth',3,'color',[.2, .6, 1])
%     hold on
%     plot(pw.full_trial.parameters.frex_2_12_Hz,data_2_plot_1{3,2}(3,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     ylim([0 .15])
%     title({['\fontsize{14}Trial ' num2str(ii)],[]})
% end
% 
% subplot(6,12,60);
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,2}(3,:,best_trials1),3),'linewidth',3,'color',[.2, .6, 1])
% hold on
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{3,2}(2,:,best_trials1),3),'linewidth',3,'color',[.7 .7 .7])
% set(gca, 'YTick', [])
% xlim([3 12])
% ylim([0 .15])
% 
% subplot(6,12,72);
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{2,2}(3,:,best_trials2),3),'linewidth',3,'color',[.2, .6, 1])
% hold on
% plot(pw.full_trial.parameters.frex_2_12_Hz, mean(data_2_plot_1{3,2}(3,:,best_trials2),3),'linewidth',3,'color',[.7 .7 .7])
% set(gca, 'YTick', [])
% xlim([3 12])
% ylim([0 .15])
% 
% 
% 
% clear ('ii',"data_2_plot_2","data_2_plot_3","data_2_plot_1")

%% Save

newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr);

saveas(gcf,name,'png')

set(gcf,'renderer','Painters')
exportgraphics(gcf,strcat(name,'.eps'),'Resolution', 300)

close all

clear('name','newStr1','path')


%% Delete bad trials

% for jj = 2:size(pw.full_trial.Pxx,1)
%     for ii = 1:size(pw.full_trial.Pxx,2)
% 
%         pw.full_trial.Pxx{jj,ii} = pw.full_trial.Pxx{jj,ii}(:,:,CSIT{ms});
% 
%         pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{jj,ii} = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{jj,ii}(:,:,CSIT{ms});
%         pw.full_trial.Pxx_baseline_norm_2_12_Hz{jj,ii}   = pw.full_trial.Pxx_baseline_norm_2_12_Hz{jj,ii}(:,:,CSIT{ms});
%         pw.full_trial.Pxx_z_norm_2_12_Hz{jj,ii}          = pw.full_trial.Pxx_z_norm_2_12_Hz{jj,ii}(:,:,CSIT{ms});
% 
%     end
% end

%% STATS - Normalization and trials average -> 3 - 5 Hz

pw.full_trial.stats = [];


% Frequency resolution. According to the fft time window
steps                                           = diff(pw.full_trial.parameters_freq_v_2_12_);

% Set freuency range
pw.full_trial.stats.parameters.frex_3_5_Hz      = 3:steps(1):5;
pw.full_trial.stats.parameters.frex_idx_3_5_Hz  = dsearchn(pw.full_trial.parameters_freq_v_2_12_,pw.full_trial.stats.parameters.frex_3_5_Hz');
pw.full_trial.stats.parameters.freq_v_3_5_      = pw.full_trial.parameters_freq_v_2_12_(pw.full_trial.stats.parameters.frex_idx_3_5_Hz);



for jj = 1:size(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz,1)
    for ii = 1:size(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz,2)

        if isempty(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{jj,ii})
            continue
        end

        % Total power normalization
        pw.full_trial.stats.Pxx_TotalPower_norm_3_5_Hz{jj,ii} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{jj,ii}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2,'omitnan'));
        
        % Baseline Normalization
        pw.full_trial.stats.Pxx_baseline_norm_3_5_Hz{jj,ii} = squeeze(mean(pw.full_trial.Pxx_baseline_norm_2_12_Hz{jj,ii}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2,'omitnan'));
        
        % Zscore Normalization
        pw.full_trial.stats.Pxx_z_norm_3_5_Hz{jj,ii} = squeeze(mean(pw.full_trial.Pxx_z_norm_2_12_Hz{jj,ii}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2,'omitnan'));
    end
end


%CS mean channels-> colunm 1 to 3 
% For now, the average of all trials will come from this normalization.
pw.full_trial.stats.Pxx_on_off_3_5_Hz{1,1} = squeeze(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,1}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2,'omitnan'));
pw.full_trial.stats.Pxx_on_off_3_5_Hz{1,2} = squeeze(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,2}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2,'omitnan'));
pw.full_trial.stats.Pxx_on_off_3_5_Hz{1,3} = squeeze(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,3}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2,'omitnan'));

pw.full_trial.stats.Pxx_on_off_3_5_Hz{2,1} = mean(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,1}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2),3,'omitnan');
pw.full_trial.stats.Pxx_on_off_3_5_Hz{2,2} = mean(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,2}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2),3,'omitnan');
pw.full_trial.stats.Pxx_on_off_3_5_Hz{2,3} = mean(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,3}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2),3,'omitnan');

%ITI mean channels -> colunm 4
%For now, I haven't calculated the values for the ITI considering all the channels—only for the average of the channels
pw.full_trial.stats.Pxx_on_off_3_5_Hz{1,4} = squeeze(mean(pw.full_trial.Pxx_on_off_2_12_Hz{2,3}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2,'omitnan'));
pw.full_trial.stats.Pxx_on_off_3_5_Hz{2,4} = mean(mean(pw.full_trial.Pxx_on_off_2_12_Hz{2,3}(:,pw.full_trial.stats.parameters.frex_idx_3_5_Hz,:),2),3,'omitnan');


% STATS - Normalization and trials average -> 7 - 9 Hz


pw.full_trial.stats.parameters.frex_7_9_Hz      = 7.8:steps(1):8.2;
pw.full_trial.stats.parameters.frex_idx_7_9_Hz  = dsearchn(pw.full_trial.parameters_freq_v_2_12_,pw.full_trial.stats.parameters.frex_7_9_Hz');
pw.full_trial.stats.parameters.freq_v_7_9_      = pw.full_trial.parameters_freq_v_2_12_(pw.full_trial.stats.parameters.frex_idx_7_9_Hz);



for jj = 1:size(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz,1)
    for ii = 1:size(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz,2)

        if isempty(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{jj,ii})
            continue
        end

        % Total power normalization
        pw.full_trial.stats.Pxx_TotalPower_norm_7_9_Hz{jj,ii} = squeeze(mean(pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{jj,ii}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2,'omitnan'));
        
        % Baseline Normalization
        pw.full_trial.stats.Pxx_baseline_norm_7_9_Hz{jj,ii} = squeeze(mean(pw.full_trial.Pxx_baseline_norm_2_12_Hz{jj,ii}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2,'omitnan'));
        
        % Zscore Normalization
        pw.full_trial.stats.Pxx_z_norm_7_9_Hz{jj,ii} = squeeze(mean(pw.full_trial.Pxx_z_norm_2_12_Hz{jj,ii}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2,'omitnan'));
    end
end



%CS mean channels-> colunm 1 to 3 
% For now, the average of all trials will come from this normalization.
pw.full_trial.stats.Pxx_on_off_7_9_Hz{1,1} = squeeze(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,1}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2,'omitnan'));
pw.full_trial.stats.Pxx_on_off_7_9_Hz{1,2} = squeeze(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,2}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2,'omitnan'));
pw.full_trial.stats.Pxx_on_off_7_9_Hz{1,3} = squeeze(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,3}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2,'omitnan'));

pw.full_trial.stats.Pxx_on_off_7_9_Hz{2,1} = mean(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,1}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2),3,'omitnan');
pw.full_trial.stats.Pxx_on_off_7_9_Hz{2,2} = mean(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,2}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2),3,'omitnan');
pw.full_trial.stats.Pxx_on_off_7_9_Hz{2,3} = mean(mean(pw.full_trial.Pxx_on_off_2_12_Hz{1,3}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2),3,'omitnan');

%ITI mean channels -> colunm 4
%For now, I haven't calculated the values for the ITI considering all the channels—only for the average of the channels
pw.full_trial.stats.Pxx_on_off_7_9_Hz{1,4} = squeeze(mean(pw.full_trial.Pxx_on_off_2_12_Hz{2,3}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2,'omitnan'));
pw.full_trial.stats.Pxx_on_off_7_9_Hz{2,4} = mean(mean(pw.full_trial.Pxx_on_off_2_12_Hz{2,3}(:,pw.full_trial.stats.parameters.frex_idx_7_9_Hz,:),2),3,'omitnan');


clear('ii','jj','steps')

%% Save

newStr = id(1:end-5);

path = '/Users/flavio/Desktop';
%path = files.FilesLoaded{1, 1}.folder;

name = strcat(path,'/',newStr,'_pw_Power');

% save data
save(name,'pw','id','-v7.3')

clear('name','newStr','path')


%% last update 20/05/2024
%  listening: Sonic Youth - Disconnection Notice
