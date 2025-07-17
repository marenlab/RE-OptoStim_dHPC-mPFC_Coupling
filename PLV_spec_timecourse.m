
disp('Filtering and Hilbert transform loop on CS-Trials')

PhaseFreqVector = 2:.5:13;
PhaseFreq_BandWidth = 1;

Data_phase = [];
Phase_diff = [];

for jj = 1:size(data_cut,1)
    for ii=1:length(PhaseFreqVector)
        Pf1 = PhaseFreqVector(ii);
        Pf2 = Pf1 + PhaseFreq_BandWidth;
        PhaseFreq = eegfilt(data_cut(jj,:),parameters.decimated_srate,Pf1,Pf2); % this is just filtering
        Data_phase(ii,:,jj) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
    end
end


Phase_diff(:,:,1) =  exp(1i*(Data_phase(:,:,3) - Data_phase(:,:,1))); % dHPC - mPFC IL
Phase_diff(:,:,2) =  exp(1i*(Data_phase(:,:,3) - Data_phase(:,:,2))); % dHPC - mPFC PL

clear('PhaseFreq_BandWidth','jj','ii','Pf1','Pf2','PhaseFreq')

%% Sliding moving average

time_bins     = 2; % seconds
time_bins_idx = time_bins * parameters.decimated_srate;

% Overlap
timeoverlap    = .95; % percentage
overlap = round((time_bins_idx)-(timeoverlap*time_bins_idx));
bins = (2:overlap:size(Phase_diff,2)-time_bins_idx);


Phase_diff_smooth = [];

for bb = 2:size(bins,2)
    Phase_diff_smooth(:,bb-1,:) = mean(Phase_diff(:,bins(bb-1):bins(bb-1) + time_bins_idx -1,:),2); % euler_differences already extracted from section above
end

PLV_spec = [];
PLV_spec = abs(Phase_diff_smooth);
time_v = linspace(data.timev_decimated(1,1),data.timev_decimated(1,end),size(PLV_spec,2));


clear('time_bins','timeoverlap','overlap','bins','bb')


%% PLot

% load python path
py_path = "/opt/anaconda3/bin/python";


title_ = {'dHPC <-> mPFC IL';'dHPC <-> mPFC PL'};

figure;%('WindowState','maximized');
set(gcf,'color','w');
sc = [1529,-366,2992,340];
set(gcf, 'Position', sc);

for ii = 1:size(PLV_spec,3)

    subplot(1,2,ii)
    contourf(time_v,PhaseFreqVector,PLV_spec(:,:,ii), 80,'linecolor','none');
    c = colorbar;
    clim([0 1])

    title(['channel ',num2str(ii+1)]);
    title(title_{ii})

    a = gca;
    a.TitleHorizontalAlignment = 'left';


    xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
    xlim([time_v(1) time_v(end)])

    clim([0.2 1])
    c = colorbar;
    c.Label.String = 'PLV';
    Py_map = getPyPlot_cMap('Blues', [], [], py_path); % chR2
    %Py_map = getPyPlot_cMap('Reds', [], [], py_path); % mCherry
    colormap(Py_map)

end

clear('py_path','title_','sc','ii','c','a','Py_map')

%% Save

clear('data')

newStr = id;
path = '/Users/flavio/Desktop';
name = strcat(path,'/',newStr,'_PLV_spect_chR2_');
%name = strcat(path,'/',newStr,'_mCherry_spect_chR2_');

save(name,'-v7.3')
saveas(gcf,name,'png')
exportgraphics(gcf,strcat(name,'.eps'))


close all



