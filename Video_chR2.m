
trial_ = 1;

t_pre = 10000;
t_pos = 10000;
stim_time = data.events{2, 1}(trial_,2) - data.events{2, 1}(trial_,1);

data_plot = data.lfp{5,1}(:,(data.events{2, 1}(trial_,1) - t_pre) : (data.events{2, 1}(trial_,2)) + t_pos);


% Choose data set
parameters.filter.longcutoff_1 = [7.5 8.5];
params.bandstop = 0;

data_plot_filtered     = [];
data_plot_phases       = [];
data_plot_phasesdiff   = [];
euler_diff             = [];

for ii = 1:3
    data_plot_filtered(ii,:) = fun_myfilters(data_plot(ii,:),parameters.decimated_srate,parameters.filter.longcutoff_1,'eegfilt',params);
    data_plot_phases(ii,:) = angle(hilbert(data_plot_filtered(ii,:)));
end


%phase_diff = data_plot_phases(3,:)-data_plot_phases(1,:);

% dHPC <--> mPFC-IL
euler_diff(1,:)             = exp(1i*(data_plot_phases(3,:)-data_plot_phases(1,:)));
data_plot_phasesdiff(1,:)   = angle(euler_diff(1,:));

% dHPC <--> mPFC-PL
euler_diff(2,:)             = exp(1i*(data_plot_phases(3,:)-data_plot_phases(2,:)));
data_plot_phasesdiff(2,:)   = angle(euler_diff(2,:));

timev = linspace(-t_pre/1000, t_pos/1000 + stim_time/1000, length(euler_diff));
center_freq = 8;

%%
time_bins     = 1; % seconds
time_bins_idx = time_bins * parameters.decimated_srate;

% Overlap
timeoverlap    = .95; % percentage
overlap = round((time_bins_idx)-(timeoverlap*time_bins_idx));


bins = (2:overlap:size(euler_diff,2)-time_bins_idx);

euler_diff_bins = [];
data_plot_phasesdiff_bins = [];
data_plot_angle_bins = [];
data_plot_PLV = [];

for ii = 1:size(euler_diff,1)
    for bb = 2:size(bins,2)


        euler_diff_bins(ii,bb-1)             = mean(euler_diff(ii,bins(bb-1):bins(bb-1) + time_bins_idx -1),2);

        data_plot_phasesdiff_bins(ii,bb-1)   = circ_mean(data_plot_phasesdiff(ii,bins(bb-1):bins(bb-1) + time_bins_idx -1),[],2);
        data_plot_PLV(ii,bb-1)               = abs(euler_diff_bins(ii,bb-1));

    end
end

data_plot_spect= [];
for ii = 1:size(data_plot,1)
    for bb = 2:size(bins,2)

        [data_plot_spect(bb-1,:,ii),hz] = pwelch(data_plot(ii,bins(bb-1):bins(bb-1) + time_bins_idx -1),1000,250,2^15,1000);

    end
end


timev_bins = linspace(-t_pre/1000, t_pos/1000 + stim_time/1000, length(euler_diff_bins));

%% Choose band F

% Frequency resolution. According to the fft time window
steps = diff(hz);

% Set freuency range
frex_      = 3:steps(1):12;
frex_idx_  = dsearchn(hz,frex_' );
freq_v     = hz(frex_idx_);

data_plot_spect_band = [];

scale_= 1/0.153;
data_plot_spect_band = normalize(scale_.*(data_plot_spect(:,frex_idx_,:)./sum(data_plot_spect(:,frex_idx_,:),2)),2,'range');
        

%%

data_plot_PLV_interp = [];
xq1 = linspace(1,length(data_plot_PLV),length(euler_diff));

for ii = 1:size(euler_diff,1)
    data_plot_PLV_interp(ii,:) = interp1(1:length(data_plot_PLV),data_plot_PLV(ii,:),xq1,'spline');
end



data_plot_spect_band_interp = [];
xq1 = linspace(1,size(data_plot_spect_band,1),length(data_plot));

for ii = 1:size(data_plot_spect_band,3)
    data_plot_spect_band_interp(:,:,ii) = interp1(1:size(data_plot_spect_band,1),data_plot_spect_band(:,:,ii),xq1,'spline');
end

%% setup figure and define plot handles


figure
set(gcf,'color','white')
set(gcf, 'Position', [1779, 266, 1854,532]);

% draw the filtered signals
subplot(6,4,[1 5])
filterplotH1 = plot(timev_bins(1),data_plot_filtered(1,1),'color',[.2, 0, .8],'linew',1);
set(gca,'xlim',[timev_bins(1) timev_bins(end)],'ylim',[min(data_plot_filtered(:)) max(data_plot_filtered(:))],'XColor', 'none')
ylabel({'mPFC IL';'\muV'})
xline(0,'--')
xline(10,'--')
xline(-5,'--')
xline(15,'--')
title('Laser on     Laser on+CS     Laser on','Color','blue')
box off

subplot(6,4,[9 13])
filterplotH2 = plot(timev_bins(1),data_plot_filtered(2,1),'color',[.2, 0, .8],'linew',1);
set(gca,'xlim',[timev_bins(1) timev_bins(end)],'ylim',[min(data_plot_filtered(:)) max(data_plot_filtered(:))],'XColor', 'none')
ylabel({'mPFC PL';'\muV'})
xline(0,'--')
xline(10,'--')
xline(-5,'--')
xline(15,'--')
box off

subplot(6,4,[17 21])
filterplotH3 = plot(timev_bins(1),data_plot_filtered(3,1),'color',[0, .6, 1],'linew',1);
set(gca,'xlim',[timev_bins(1) timev_bins(end)],'ylim',[min(data_plot_filtered(:)) max(data_plot_filtered(:))])
xlabel('Time (s)')
ylabel({'dHPC';'\muV'})
xline(0,'--')
xline(10,'--')
xline(-5,'--')
xline(15,'--')
box off



subplot(6,4,[2 6])
spectH1 = plot(freq_v,zeros(length(freq_v)),'color',[.2, 0, .8],'linew',2);
set(gca,'xlim',[3 12],'ylim',[0 1],'XColor', 'none');
ylabel({'Amplitude (a.u.)';[]})
box off

subplot(6,4,[10 14])
spectH2 = plot(freq_v,zeros(length(freq_v)),'color',[.2, 0, .8],'linew',2);
set(gca,'xlim',[3 12],'ylim',[0 1],'XColor', 'none');
ylabel({'Amplitude (a.u.)';[]})
box off

subplot(6,4,[18 22])
spectH3 = plot(freq_v,zeros(length(freq_v)),'color',[0, .6, 1],'linew',2);
set(gca,'xlim',[3 12],'ylim',[0 1]);
xlabel('Frequency (Hz)')
ylabel({'Amplitude (a.u.)';[]})
box off



subplot(6,4,[3 11])
phaseanglesH1 = plot(timev(1),data_plot_phases(1,1),'color',[.2, 0, .8],'linew',2);
hold on
phaseanglesH2 = plot(timev(1),data_plot_phases(3,1),'color',[0, .6, 1],'linew',2);
set(gca,'xlim',[timev(1) timev(end)],'ylim',[min(data_plot_phases(:)) max(data_plot_phases(:))],'XColor', 'none')
ylabel({'hHPC - mPFC IL';'Degrees (radians)';[]})
ax1 = gca;
box off

subplot(6,4,[15 23]);
phaseanglesH3 = plot(timev(1),data_plot_phases(2,1),'color',[.2,  0, .8],'linew',2);
hold on
phaseanglesH4 = plot(timev(1),data_plot_phases(3,1),'color',[0, .6, 1],'linew',2);
set(gca,'xlim',[timev(1) timev(end)],'ylim',[min(data_plot_phases(:)) max(data_plot_phases(:))])
ylabel({'hHPC - mPFC PL';'Degrees (radians)';[]})
xlabel('Time (s)')
ax2 = gca;
box off

s1 = subplot(6,4,[4 12]);
s1.Position =  [0.7484    0.5374    0.9*0.1566    0.9*0.3876];
polarAngleDiffH1 = polarplot([zeros(1,1) data_plot_phasesdiff_bins(1,1)]',repmat([0 1],1,1)');

% hold on
% polarAngleDiffH1 = polarplot([zeros(1,1) data_plot_phasesdiff_bins(1,1)]',repmat([0 1],1,1)','color',[0, .6, 1],'linew', 3);
% legend('phase difference','phase lock value')

s2 = subplot(6,4,[16 24]);
s2.Position =  [0.7484    0.1100    0.9*0.1566    0.9*0.3876];

polarAngleDiffH2 = polarplot([zeros(1,1) data_plot_phasesdiff_bins(2,1)]',repmat([0 1],1,1)');
% hold on
% polarAngleDiffH2 = polarplot([zeros(1,1) data_plot_phasesdiff_bins(1,1)]',repmat([0 1],1,1)','color',[0, .6, 1],'linew', 3);
% legend('phase difference','phase lock value')

    annotation(gcf,'textbox',...
        [0.871010787486511 0.417293233082706 0.0745124056094955 0.0451127819548869],...
        'String','dHPC <---> mPFC- PL',...
        'FontSize',12,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);

    % Create textbox
    annotation(gcf,'textbox',...
        [0.872628910463858 0.857142857142856 0.0707367853290192 0.045112781954887],...
        'String','dHPC <---> mPFC- IL',...
        'FontSize',12,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);

    % Create textbox
    annotation(gcf,'textbox',...
        [0.923869471413153 0.060150375939847 0.0745124056094987 0.0451127819548869],...
        'String',{'Phase difference','Phase lock value'},...
        'FontSize',12,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);

    % Create line
    annotation(gcf,'line',[0.899137001078749 0.920711974110032],...
        [0.0854661654135338 0.0845864661654135]);

    % Create line
    annotation(gcf,'line',[0.89967637540453 0.921251348435813],...
        [0.0572706766917293 0.056390977443609],...
        'Color',[.2, 0, .8],...
        'LineWidth',3);

%% now update plots at each timestep


mov = VideoWriter('chR2','Uncompressed AVI');

open(mov)


for ti=2:10:length(data_plot)
    
    % update filtered signals
    set(filterplotH1,'XData',timev(1:ti),'YData',data_plot_filtered(1,1:ti))
    set(filterplotH2,'XData',timev(1:ti),'YData',data_plot_filtered(2,1:ti))
    set(filterplotH3,'XData',timev(1:ti),'YData',data_plot_filtered(3,1:ti))

    % % update fft
    set(spectH1,'YData',mean(data_plot_spect_band_interp(ti:ti+10,:,1),1))
    set(spectH2,'YData',mean(data_plot_spect_band_interp(ti:ti+10,:,2),1))
    set(spectH3,'YData',mean(data_plot_spect_band_interp(ti:ti+10,:,3),1))


   
    % update cartesian plot of phase angles
    set(phaseanglesH1,'XData',timev(1:ti),'YData',data_plot_phases(1,1:ti))
    set(phaseanglesH2,'XData',timev(1:ti),'YData',data_plot_phases(3,1:ti))
    ax1.XLim = ([timev(ti)-1 timev(ti)]);   

    set(phaseanglesH3,'XData',timev(1:ti),'YData',data_plot_phases(2,1:ti))
    set(phaseanglesH4,'XData',timev(1:ti),'YData',data_plot_phases(3,1:ti))    
    ax2.XLim = ([timev(ti)-1 timev(ti)]);


    s1 = subplot(6,4,[4 12]);
    s1.Position =  [0.7484    0.5374    0.9*0.1566    0.9*0.3876];
    cla
    polarplot([zeros(1,ti) data_plot_phasesdiff(1,1:ti)]',repmat([0 1],1,ti)','color',[.6 .6 .6 .4]);
    hold on
    polarplot([0 circ_mean(data_plot_phasesdiff(1,ti:ti+1),[],2)],[0 mean(data_plot_PLV_interp(1,ti:ti+1)-.2)],'color',[.2, 0, .8],'linew', 5);

    s2 = subplot(6,4,[16 24]);
    s2.Position =  [0.7484    0.1100    0.9*0.1566    0.9*0.3876];
    cla
    polarplot([zeros(1,ti) data_plot_phasesdiff(2,1:ti)]',repmat([0 1],1,ti)','color',[.6 .6 .6 .4]);
    hold on
    polarplot([0 circ_mean(data_plot_phasesdiff(2,ti:ti+1),[],2)],[0 mean(data_plot_PLV_interp(2,ti:ti+1)-.2)],'color',[.2, 0, .8],'linew', 5);



    drawnow

    F = getframe(gcf);
    writeVideo(mov,F);

end

close(mov)


%%