
%% Continuous 1-D wavelet transform

% by Flavio Mourao.
% email: mourao.fg@illinois.edu
% Maren Lab -  Beckman Institute for Advanced Science and Technology
% TUniversity of Illinois Urbana-Champaign

% Started in:  05/2025
% Last update: 06/2025

%%
% First extract the data with: Extracting_LFPs_and_events.m
% ... then organize the data with the script: Pre_processing.m

%%

fprintf('\n Continuous 1-D wavelet transform... \n');

%% Normalize data

data_z = normalize(data.lfp{5,1});
% data_z = normalize(data.lfp{5,1}(:,data.events{2, 1}(1,1)-20000 : data.events{2, 1}(5,2)+20000),2);


%% CWT Matlab function

%lines: frequencies / columns: time / third dimension: channels

wavelt = [];

% Continuous wavelet transform filter bank
%fb = cwtfilterbank('SignalLength',sigLength,'SamplingFrequency',Fs,'FrequencyLimits',[2 12],'Wavelet','morse','TimeBandwidth',20);

% freqz(fb) % Check gaussians

waveletParameters = [2,80];

for ii = 1:size(data_z,1)
    if ii ==1
       [wavelt.data{1,1}(:,:,ii),wavelt.freq] = cwt(data_z(ii,:),parameters.original_srate, WaveletParameters = waveletParameters,FrequencyLimits = [2 12]);
       %[wavelt.data{1,1}(:,:,ii),wavelt.freq] = cwt(data_z(ii,:),parameters.original_srate,'amor',FrequencyLimits = [2 12]);

    else
        wavelt.data{1,1}(:,:,ii) = cwt(data_z(ii,:),parameters.original_srate, WaveletParameters = waveletParameters,FrequencyLimits = [2 12]);
        %wavelt.data{1,1}(:,:,ii) = cwt(data_z(ii,:),parameters.original_srate,'amor',FrequencyLimits = [2 12]);
    end
end


wavelt.time = (0:size(data_z,2)-1)/parameters.original_srate;

clear('ii')

%% checking other channels

% PL_  = [2];   % mPFC IL
% IL_  = [15];  % mPFC PL
% dHPC = [17];  % dHPC
% 
% wavelt.data{1, 2}(:,:,1) = mean(abs(wavelt.data{1, 1}(:,:,IL_)),3);   % mPFC-IL
% wavelt.data{1, 2}(:,:,2) = mean(abs(wavelt.data{1, 1}(:,:,PL_)),3);   % mPFC-PL
% wavelt.data{1, 2}(:,:,3) = mean(abs(wavelt.data{1, 1}(:,:,dHPC)),3);  % dHPC
% 
% clear('PL_','IL_','dHPC')

%% Plot to check full session. Channels per substrate 

data_to_plot = wavelt.data{1,1};

% load python path
py_path = "/opt/anaconda3/bin/python";


title_ = {'mPFC IL';'mPFC PL';'dHPC'};

figure;%('WindowState','maximized');
set(gcf,'color','w');
sc = [1529,-366,2992,340];
set(gcf, 'Position', sc);

for ii = 1:size(data_to_plot,3)

    subplot(1,size(data_to_plot,3),ii)
    contourf(wavelt.time(1:250:end),wavelt.freq(:,1),zscore(abs(data_to_plot(:,1:250:end,ii)),[],1), 80,'linecolor','none'); % 250 ms step / z scored Freq domain

    title(['channel ',num2str(ii+1)]);
    title(title_{ii})

    a = gca;
    a.TitleHorizontalAlignment = 'left';

    xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14)
    xlim([wavelt.time(1) wavelt.time(end)])

    clim([0 3])
    c = colorbar;
    c.Label.String = 'Normalized Power (A.U.)';
    %Py_map = getPyPlot_cMap('Blues', [], [], py_path); % chR2
    Py_map = getPyPlot_cMap('Reds', [], [], py_path);   % mCherry
    colormap(Py_map)

end

clear ('data_to_plot','ii','a','sc','c','py_path','title_')

%% Save

newStr1 = id(1:end-12);

% %path = files.FilesLoaded{1, 1}.folder;
path = '/Users/flavio/Desktop';
name = strcat(path,'/',newStr1,'_wavelt');
saveas(gcf,name,'png')

% set(gcf,'renderer','Painters')
% exportgraphics(gcf,strcat(name,'.eps'),'Resolution', 300)

save(name,'wavelt','id','-v7.3')

close all

clear('name','newStr1','path')
%% last update 05/26/2025
%  listening: American Football : for sure

