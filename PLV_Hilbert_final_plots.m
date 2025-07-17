% Phase Coherence (Hilbert Transform (Instantaneous phase)) final plot

% by Flavio Mourao.
% email: mourao.fg@tamu.edu
% Maren Lab - Department of Psychological and Brain Sciences
% Texas A&M University

% Started in:  05/2024
% Last update: 07/2025

%% Prepare data and Load data

%% Extinction/Retrieval

for ii = 1:size(hilb.delta_phase,1)
    for    jj = 1:size(hilb.delta_phase,2)

        data_phase{ii,jj}(:,:,:,ms)       = hilb.delta_phase_5_trial{ii,jj};
        data_phase_mean{ii,jj}(:,:,:,ms)  = hilb.delta_phase_mean_5_trial{ii,jj};
        data_PLV_mean{ii,jj}(:,:,:,ms)    = hilb.PLV_mean_5_trial{ii,jj};



    end
end

clear('ii','jj')


% data_phase{p,2} = data_phase{p,2}(:,:,:,[1 2 4]);

%%
if ms == length(files.FilesLoaded{1, 1})

    % Mean and SEM
    data_2_plot_data_phase      = cellfun(@(x)circ_mean(x,[],4),data_phase,'UniformOutput',false);
    data_2_plot_data_phase_mean = cellfun(@(x)circ_mean(x,[],4),data_phase_mean,'UniformOutput',false);
    data_2_plot_data_PLV_mean   = cellfun(@(x)median(x,4),data_PLV_mean,'UniformOutput',false);

    %% Polar PLots. Baseline and Mean Trials. 4 Hz and 8 Hz

p = 2; % 2 - CS-trials / 3 - ITI

% figure
hold on
color_ = [.2, .6, 1, .1]; % CHR2
color_m = [0, 0, .6];

% color_ = [.8, 0, 0, .1]; % mCherry
% color_m = [0.6350, 0.0780, 0.1840];

% Define data to plot
phase = data_2_plot_data_phase;

phase_mean = data_2_plot_data_phase_mean  ;
PLV_mean = data_2_plot_data_PLV_mean ;

sgtitle_ = {'4 Hz Bandpass filter';'8 Hz Bandpass filter'};

% Baseline
for ii = 1:size(phase,2)

    figure(ii);
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','w');

    sgtitle({sgtitle_{ii};''},'FontWeight','bold','FontSize',16);

    % mPFC PL <--> mPFC IL
    subplot(3,10,1)
    polarplot([zeros(1,size(phase{1,ii}(1,1:240:end),2)) phase{1,ii}(1,1:240:end)]',repmat([0 1],1,size(phase{1,ii}(1,1:240:end),2))','color',color_);
    hold on
    I=polarplot([0 phase_mean{1,ii}(1,1)],[0 PLV_mean{1,ii}(1,1)]);
    set(I,'linewidth',6,'color',color_m)
    thetaticks([0 90 180 270])
    title({['\fontsize{16} Baseline'];[];['\fontsize{14} PLV = ', num2str(PLV_mean{1,ii}(1,1))];[]})

    % mPFC IL <--> dHPC
    subplot(3,10,11)
    polarplot([zeros(1,size(phase{1,ii}(2,1:240:end),2)) phase{1,ii}(2,1:240:end)]',repmat([0 1],1,size(phase{1,ii}(2,1:240:end),2))','color',color_);
    hold on
    I=polarplot([0 phase_mean{1,ii}(2,1)],[0 PLV_mean{1,ii}(2,1)]);
    set(I,'linewidth',6,'color',color_m)
    thetaticks([0 90 180 270])
    title({['\fontsize{14} PLV = ', num2str(PLV_mean{1,ii}(2,1))];[]})

    % mPFC PL <--> dHPC
    subplot(3,10,21)
    polarplot([zeros(1,size(phase{1,ii}(3,1:240:end),2)) phase{1,ii}(3,1:240:end)]',repmat([0 1],1,size(phase{1,ii}(3,1:240:end),2))','color',color_);
    hold on
    I=polarplot([0 phase_mean{1,ii}(3,1)],[0 PLV_mean{1,ii}(3,1)]);
    set(I,'linewidth',6,'color',color_m)
    thetaticks([0 90 180 270])
    title({['\fontsize{14} PLV = ', num2str(PLV_mean{1,ii}(3,1))];[]})


    annotation(figure(ii),'textbox',...
        [0.0172619047619056 0.807610993657511 0.101190476190476 0.0243699788583516],...
        'String','mPFC PL <--> mPFC IL',...
        'FontWeight','bold',...
        'FontSize',14,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);


    annotation(figure(ii),'textbox',...
        [0.022619047619049 0.524312896405924 0.101190476190477 0.0243699788583518],...
        'String','mPFC PL <--> dHPC',...
        'FontWeight','bold',...
        'FontSize',14,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);

    annotation(figure(ii),'textbox',...
        [0.0267857142857152 0.250528541226218 0.101190476190476 0.0243699788583517],...
        'String','mPFC IL <--> dHPC',...
        'FontWeight','bold',...
        'FontSize',14,...
        'FontName','Helvetica Neue',...
        'FitBoxToText','off',...
        'EdgeColor',[1 1 1]);



    % annotation(figure(ii),'textbox',...
    % [0.0148809523809532 0.957716701902754 0.101190476190476 0.0243699788583516],...
    % 'String',id(1:end-11),...
    % 'FontWeight','bold',...
    % 'FontSize',18,...
    % 'FontName','Helvetica Neue',...
    % 'FitBoxToText','off',...
    % 'EdgeColor',[1 1 1]);

end

% Mean 8 Hz
for ii = 1:size(phase,2)
    figure(ii)
    for tt = 1:size(phase{p,ii},3)

        % mPFC PL <--> mPFC IL
        subplot(3,10,tt+1)
        polarplot([zeros(1,size(phase{p,ii}(1,1:240:end,tt),2)) phase{p,ii}(1,1:240:end,tt)]',repmat([0 1],1,size(phase{p,ii}(1,1:240:end,tt),2))','linewidth',.5,'color',color_);
        hold on
        I=polarplot([0 phase_mean{p,ii}(1,tt)],[0 PLV_mean{p,ii}(1,tt)]);
        set(I,'linewidth',6,'color',color_m)
        thetaticks([0 90 180 270])
        title({['\fontsize{16} 5 min Block ' num2str(tt)];[];['\fontsize{14} PLV = ' num2str(PLV_mean{4,ii}(1,tt))];[]},'color',color_)

        % mPFC PL <--> dHPC
        subplot(3,10,tt+11)
        polarplot([zeros(1,size(phase{p,ii}(2,1:240:end,tt),2)) phase{p,ii}(2,1:240:end,tt)]',repmat([0 1],1,size(phase{p,ii}(2,1:240:end,tt),2))','linewidth',.5,'color',color_);
        hold on
        I=polarplot([0 phase_mean{p,ii}(2,tt)],[0 PLV_mean{p,ii}(2,tt)]);
        set(I,'linewidth',6,'color',color_m)
        thetaticks([0 90 180 270])
        title({['\fontsize{14} PLV = ', num2str(PLV_mean{p,ii}(2,tt))];[]},'color',color_)

        % mPFC IL <--> dHPC
        subplot(3,10,tt+21)
        polarplot([zeros(1,size(phase{p,ii}(3,1:240:end,tt),2)) phase{p,ii}(3,1:240:end,tt)]',repmat([0 1],1,size(phase{p,ii}(3,1:240:end,tt),2))','linewidth',.5,'color',color_);
        hold on
        I=polarplot([0 phase_mean{p,ii}(3,tt)],[0 PLV_mean{p,ii}(3,tt)]);
        set(I,'linewidth',6,'color',color_m)
        thetaticks([0 90 180 270])
        title({['\fontsize{14} PLV = ', num2str(PLV_mean{p,ii}(3,tt))];[]},'color',color_)

    end
end
end


%clear('ii','I','data_phase','data_phase_mean','data_PLV_mean')



%% Save Figures
%
% newStr = id(1:end-8);
%
% path = '/Users/flavio/Desktop';
% %path = files.FilesLoaded{1, 1}.folder;
%
% name = strcat(path,'/',newStr,'_PLV_extinction');
%
% %get all figures
% figlist = findobj(allchild(0), 'flat', 'Type', 'figure');
%
% name_figs = {'_3_5Hz','_6_8Hz','_8_10Hz'};
%
% % Loop through figure
% for ii = 1:numel(figlist)
%     name_loop = strcat(name,name_figs{ii});
%     FigHandle = figlist(ii);
%     saveas(FigHandle,name_loop,'png')
% end
%
% close all
% clear('name','newStr','path')
%
