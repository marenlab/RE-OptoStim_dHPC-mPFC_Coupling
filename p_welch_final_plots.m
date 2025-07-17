% Welch`s final plot

% by Flavio Mourao.
% email: mourao.fg@illinois.edu
% Maren Lab -  Beckman Institute for Advanced Science and Technology
% TUniversity of Illinois Urbana-Champaign

% Started in:  05/2023
% Last update: 05/2023


%% Load data Extinction

% Habituation
% data_2_plot_all{1,ms} = pw.full_trial.Pxx_on_off_2_12_Hz{1,2}; % Stim
% data_2_plot_all{2,ms} = pw.full_trial.Pxx_on_off_2_12_Hz{2,2}; % ITI

% Habituation
data_2_plot_all{1,ms} = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{2,2}; % Stim
data_2_plot_all{2,ms} = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{3,2}; % ITI
data_2_plot_all{3,ms} = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz{1,2}; % CS-Trial Blk 


% Extinction
% data_2_plot_all{1,ms} = pw.full_trial.Pxx_on_off_2_12_Hz{1,3}; % CS+Stim
% data_2_plot_all{2,ms} = pw.full_trial.Pxx_on_off_2_12_Hz{2,3}; % ITI

% Extinction
% data_2_plot_all{1,ms} = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz  {2,3}; % CS+Stim
% data_2_plot_all{2,ms} = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz  {3,3}; % ITI
% data_2_plot_all{3,ms} = pw.full_trial.Pxx_TotalPower_norm_2_12_Hz  {1,2}; % CS-Trial Blk 


% Extinction
% data_2_plot_all{1,ms} = pw.full_trial.Pxx_z_norm_2_12_Hz{2,3}; % CS+Stim
% data_2_plot_all{2,ms} = pw.full_trial.Pxx_z_norm_2_12_Hz{3,3}; % ITI



if ms == length(files.FilesLoaded{1, 1})

    %%
    
    % for ii = 1:size(data_2_plot_all,1)
    %     for jj = 1:size(data_2_plot_all,2)
    %         data_2_plot_all{ii,jj} = normalize(data_2_plot_all{ii,jj},2,'range');
    % 
    %     end
    % end

    %% Frequency vector
    freq_v = pw.full_trial.parameters.frex_2_12_Hz;

    data_2_plot_meanTrials_SEM = [];

    for ii = 1:size(data_2_plot_all,1)

        % Concatenate animals
        data_2_plot_meanTrials_SEM{ii,1} = cat(4,(data_2_plot_all{ii,:}));

        % Mean and SEM

        data_2_plot_meanTrials_SEM{ii,2}  = median(data_2_plot_meanTrials_SEM{ii,1},4,'omitnan');
        data_2_plot_meanTrials_SEM{ii,3}  = std(data_2_plot_meanTrials_SEM{ii,1},[],4,'omitnan')./size(data_2_plot_meanTrials_SEM{ii,1},4);

    end


    data_2_plot_mean_SEM = [];

    for ii = 1:size(data_2_plot_all,1)

        % Concatenate animals
        data_2_plot_mean_SEM{ii,1} = squeeze(mean(cat(4,(data_2_plot_all{ii,:})),3,'omitnan'));

        % Mean and SEM

        data_2_plot_mean_SEM{ii,2}  = median(data_2_plot_mean_SEM{ii,1},3,'omitnan');
        data_2_plot_mean_SEM{ii,3}  = std(data_2_plot_mean_SEM{ii,1},[],3,'omitnan')./size(data_2_plot_mean_SEM{ii,1},3);

    end

    clear('ii')

%% Habituation Plot 1
% 
% data_2_plot_1      = data_2_plot_meanTrials_SEM;
% data_2_plot_2      = data_2_plot_mean_SEM;
% 
% %color_ = [.2, .6, 1]; % CHR2
% color_ = [.8, 0, 0]; % mCherry
% 
% 
% figure;
% %set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf, 'Position', [1683, -377, 1249,1068]);
% set(gcf,'color','w');
% axis('square')
% sgtitle({'\fontsize{18} \bf Habituation',[]})
% %sgtitle('\fontsize{18} \bf Retrieval')
% 
% 
% 
% % mPFC IL
% 
% for ii = 1:5
% 
%     subplot(3,6,ii);
%     boundedline(freq_v,data_2_plot_1{1,2}(1,:,ii),data_2_plot_1{1,3}(1,:,ii),'linewidth',3,'color',color_)
%     hold on 
%     boundedline(freq_v,data_2_plot_1{2,2}(1,:,ii),data_2_plot_1{2,3}(1,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     %set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     %ylim([0 0.12])
%     title({['\fontsize{14}Trial ' num2str(ii)],[]})
% end
% 
% subplot(3,6,6);
% boundedline(freq_v, data_2_plot_2{1,2}(1,:),data_2_plot_2{1,3}(1,:),'linewidth',3,'color',color_)
% hold on
% boundedline(freq_v, data_2_plot_2{2,2}(1,:),data_2_plot_2{2,3}(1,:),'linewidth',3,'color',[.7 .7 .7])
% 
% %set(gca, 'YTick', [])
% xlim([3 12])
% %ylim([0 0.12])
% title({'\fontsize{14}Averaged Trials',[]})
% 
% 
% 
% % mPFC PL
% for ii = 1:5
% 
%     subplot(3,6,ii+6);
%     boundedline(freq_v,data_2_plot_1{1,2}(1,:,ii),data_2_plot_1{1,3}(2,:,ii),'linewidth',3,'color',color_)
%     hold on 
%     boundedline(freq_v,data_2_plot_1{2,2}(1,:,ii),data_2_plot_1{2,3}(2,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     %set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     %ylim([0 0.12])
%     title({['\fontsize{14}Trial ' num2str(ii)],[]})
% end
% 
% subplot(3,6,12);
% boundedline(freq_v, data_2_plot_2{1,2}(2,:),data_2_plot_2{1,3}(2,:),'linewidth',3,'color',color_)
% hold on
% boundedline(freq_v, data_2_plot_2{2,2}(2,:),data_2_plot_2{2,3}(2,:),'linewidth',3,'color',[.7 .7 .7])
% %set(gca, 'YTick', [])
% xlim([3 12])
% %ylim([0 0.12])
% title({'\fontsize{14}Averaged Trials',[]})
% 
% 
% 
% % dHPC
% for ii = 1:5
% 
%     subplot(3,6,ii+12);
%     boundedline(freq_v,data_2_plot_1{1,2}(3,:,ii),data_2_plot_1{1,3}(3,:,ii),'linewidth',3,'color',color_)
%     hold on 
%     boundedline(freq_v,data_2_plot_1{2,2}(3,:,ii),data_2_plot_1{2,3}(3,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     %set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     %ylim([0 0.12])
%     title({['\fontsize{14}Trial ' num2str(ii)],[]})
% end
% 
% subplot(3,6,18);
% boundedline(freq_v, data_2_plot_2{1,2}(3,:),data_2_plot_2{1,3}(3,:),'linewidth',3,'color',color_)
% hold on
% boundedline(freq_v, data_2_plot_2{2,2}(3,:),data_2_plot_2{2,3}(3,:),'linewidth',3,'color',[.7 .7 .7])
% %set(gca, 'YTick', [])
% xlim([3 12])
% %ylim([0 0.12])
% title({'\fontsize{14}Averaged Trials',[]})
% 
% 
% clear ('ii',"data_2_plot_2","data_2_plot_1")
% 
% %% Habituation Plot 2
% 
% data_2_plot_1      = data_2_plot_meanTrials_SEM;
% data_2_plot_2      = data_2_plot_mean_SEM;
% 
% color_ = [.2, .6, 1]; % CHR2
% %color_ = [.8, 0, 0]; % mCherry
% 
% 
% % figure;
% %set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf, 'Position', [1683, -377, 1744,868]);
% set(gcf,'color','w');
% %axis('square')
% %sgtitle({'\fontsize{18} \bf Habituation',[]})
% %sgtitle({'\fontsize{18} \bf Extinction',[]})
% 
% sgtitle('\fontsize{18} \bf Retrieval')
% 
% hold on
% 
% % mPFC IL
% 
% for ii = 1:5
% 
%     subplot(3,6,ii);
% 
%     for aa = 1:4
%         plot(freq_v,data_2_plot_1{1,1}(1,:,ii,aa),'linewidth',1,'color',[.2, .6, 1, .5])    
%         hold on
%     end
% 
%     plot(freq_v,data_2_plot_1{1,2}(1,:,ii),'linewidth',3,'color',color_)
% 
% 
%     %set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     %%ylim([0 0.1])
%     title({['\fontsize{14}CS-Block ' num2str(ii)],[]})
% end
% 
% subplot(3,6,6);
% 
% for aa = 1:4
%     plot(freq_v, data_2_plot_2{1,1}(1,:,aa),'linewidth',1,'color',[.2, .6, 1, .5])
%     hold on
% end
% 
%     plot(freq_v, data_2_plot_2{1,2}(1,:),'linewidth',3,'color',color_)
% 
% %set(gca, 'YTick', [])
% xlim([3 12])
% %%ylim([0 0.1])
% title({'\fontsize{14}Averaged Trials',[]})
% 
% 
% 
% % mPFC PL
% for ii = 1:5
% 
%     subplot(3,6,ii+6);
% 
%     for aa = 1:4
%         plot(freq_v,data_2_plot_1{1,1}(2,:,ii,aa),'linewidth',1,'color',[.2, .6, 1, .5])    
%         hold on
%     end
% 
%     plot(freq_v,data_2_plot_1{1,2}(2,:,ii),'linewidth',3,'color',color_)
% 
% 
%     %set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     %%ylim([0 0.1])
%     title({['\fontsize{14}CS-Block ' num2str(ii)],[]})
% end
% 
% subplot(3,6,12);
% 
% for aa = 1:4
%     plot(freq_v, data_2_plot_2{1,1}(2,:,aa),'linewidth',1,'color',[.2, .6, 1, .5])
%     hold on
% end
% 
%     plot(freq_v, data_2_plot_2{1,2}(2,:),'linewidth',3,'color',color_)
% 
% %set(gca, 'YTick', [])
% xlim([3 12])
% %%ylim([0 0.1])
% title({'\fontsize{14}Averaged Trials',[]})
% 
% 
% 
% % dHPC
% for ii = 1:5
% 
%     subplot(3,6,ii+12);
% 
%     for aa = 1:4
%         plot(freq_v,data_2_plot_1{1,1}(3,:,ii,aa),'linewidth',1,'color',[.2, .6, 1, .5])    
%         hold on
%     end
% 
%     plot(freq_v,data_2_plot_1{1,2}(3,:,ii),'linewidth',3,'color',color_)
% 
% 
%     %set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     %%ylim([0 0.1])
%     title({['\fontsize{14}CS-Block ' num2str(ii)],[]})
% end
% 
% subplot(3,6,18);
% 
% for aa = 1:4
%     plot(freq_v, data_2_plot_2{1,1}(3,:,aa),'linewidth',1,'color',[.2, .6, 1, .5])
%     hold on
% end
% 
%     plot(freq_v, data_2_plot_2{1,2}(3,:),'linewidth',3,'color',color_)
% 
% %set(gca, 'YTick', [])
% xlim([3 12])
% %%ylim([0 0.1])
% title({'\fontsize{14}Averaged Trials',[]})
% 
% 
% clear ('ii',"data_2_plot_2","data_2_plot_1")

%% Extinction Plot 1

% data_2_plot_1      = data_2_plot_meanTrials_SEM;
% data_2_plot_2      = data_2_plot_mean_SEM;
% 
% %color_ = [.2, .6, 1]; % CHR2
% color_ = [.2, .6, 1]; % mCherry
% 
% 
% figure;
% %set(gcf, 'Position', get(0, 'Screensize'));
% set(gcf, 'Position', [1683, -377, 1744,868]);
% set(gcf,'color','w');
% axis('square')
% %sgtitle({'\fontsize{18} \bf Habituation',[]})
% sgtitle({'\fontsize{18} \bf Extinction',[]})
% 
% %sgtitle('\fontsize{18} \bf Retrieval')
% 
% hold on
% 
% % mPFC IL
% 
% for ii = 1:9
% 
%     subplot(3,10,ii);
%     boundedline(freq_v,data_2_plot_1{1,2}(1,:,ii),data_2_plot_1{1,3}(1,:,ii),'linewidth',3,'color',color_)
%     hold on 
%     boundedline(freq_v,data_2_plot_1{2,2}(1,:,ii),data_2_plot_1{2,3}(1,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     %set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     %ylim([0 0.12])
%     title({['\fontsize{14}CS-Block ' num2str(ii)],[]})
% end
% 
% subplot(3,10,10);
% boundedline(freq_v, data_2_plot_2{1,2}(1,:),data_2_plot_2{1,3}(1,:),'linewidth',3,'color',color_)
% hold on
% boundedline(freq_v, data_2_plot_2{2,2}(1,:),data_2_plot_2{2,3}(1,:),'linewidth',3,'color',[.7 .7 .7])
% 
% %set(gca, 'YTick', [])
% xlim([3 12])
% %%ylim([0 0.1])
% title({'\fontsize{14}Averaged Trials',[]})
% 
% 
% 
% % mPFC PL
% for ii = 1:9
% 
%     subplot(3,10,ii+10);
%     boundedline(freq_v,data_2_plot_1{1,2}(1,:,ii),data_2_plot_1{1,3}(2,:,ii),'linewidth',3,'color',color_)
%     hold on 
%     boundedline(freq_v,data_2_plot_1{2,2}(1,:,ii),data_2_plot_1{2,3}(2,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     %set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     %ylim([-0.05 3.5])
%     title({['\fontsize{14}CS-Block ' num2str(ii)],[]})
% end
% 
% subplot(3,10,20);
% boundedline(freq_v, data_2_plot_2{1,2}(2,:),data_2_plot_2{1,3}(2,:),'linewidth',3,'color',color_)
% hold on
% boundedline(freq_v, data_2_plot_2{2,2}(2,:),data_2_plot_2{2,3}(2,:),'linewidth',3,'color',[.7 .7 .7])
% %set(gca, 'YTick', [])
% xlim([3 12])
% %%ylim([0 0.1])
% title({'\fontsize{14}Averaged Trials',[]})
% 
% 
% 
% % dHPC
% for ii = 1:9
% 
%     subplot(3,10,ii+20);
%     boundedline(freq_v,data_2_plot_1{1,2}(3,:,ii),data_2_plot_1{1,3}(3,:,ii),'linewidth',3,'color',color_)
%     hold on 
%     boundedline(freq_v,data_2_plot_1{2,2}(3,:,ii),data_2_plot_1{2,3}(3,:,ii),'linewidth',3,'color',[.7 .7 .7])
%     %set(gca, 'YTick', [],'YColor','none')
%     box off
%     xlim([3 12])
%     %ylim([0 0.12])
%     title({['\fontsize{14}CS-Block ' num2str(ii)],[]})
% end
% 
% subplot(3,10,30);
% boundedline(freq_v, data_2_plot_2{1,2}(3,:),data_2_plot_2{1,3}(3,:),'linewidth',3,'color',color_)
% hold on
% boundedline(freq_v, data_2_plot_2{2,2}(3,:),data_2_plot_2{2,3}(3,:),'linewidth',3,'color',[.7 .7 .7])
% %set(gca, 'YTick', [])
% xlim([3 12])
% %%ylim([0 0.1])
% title({'\fontsize{14}Averaged Trials',[]})
% 
% 
% clear ('ii',"data_2_plot_2","data_2_plot_1")


%% Extinction Plot 2

data_2_plot_1      = data_2_plot_meanTrials_SEM;
data_2_plot_2      = data_2_plot_mean_SEM;

% color_ = [.2, .6, 1]; % CHR2
% color_1 = [.2, .6, 1, .5];

color_ = [.8, 0, 0]; % mCherry
color_1 = [.8, 0, 0, .5];
% 
p = 1;
t = 5;

figure;
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf, 'Position', [1683, -377, 1744,868]);
set(gcf,'color','w');
% axis('square')
% %sgtitle({'\fontsize{18} \bf Habituation',[]})
% sgtitle({'\fontsize{18} \bf Extinction',[]})
% 
% %sgtitle('\fontsize{18} \bf Retrieval')

hold on

% mPFC IL

for ii = 1:t

    subplot(3,10,ii);

    for aa = 1:4
        plot(freq_v,data_2_plot_1{p,1}(1,:,ii,aa),'linewidth',1,'color',[color_1])    
        hold on
    end

    plot(freq_v,data_2_plot_1{p,2}(1,:,ii),'linewidth',3,'color',color_)


    %set(gca, 'YTick', [],'YColor','none')
    box off
    xlim([3 12])
    %ylim([-0.05 3.5])
    title({['\fontsize{14}CS-Trial Blk  ' num2str(ii)],[]})
end

subplot(3,10,10);

for aa = 1:4
    plot(freq_v, data_2_plot_2{p,1}(1,:,aa),'linewidth',1,'color',[color_1])
    hold on
end

    plot(freq_v, data_2_plot_2{p,2}(1,:),'linewidth',3,'color',color_)

%set(gca, 'YTick', [])
xlim([3 12])
%ylim([0 0.1])
title({'\fontsize{14}Averaged Trials',[]})



% mPFC PL
for ii = 1:t

    subplot(3,10,ii+10);

    for aa = 1:4
        plot(freq_v,data_2_plot_1{p,1}(2,:,ii,aa),'linewidth',1,'color',[color_1])    
        hold on
    end

    plot(freq_v,data_2_plot_1{p,2}(2,:,ii),'linewidth',3,'color',color_)


    %set(gca, 'YTick', [],'YColor','none')
    box off
    xlim([3 12])
    %ylim([-0.05 3.5])
    title({['\fontsize{14}CS-Trial Blk  ' num2str(ii)],[]})
end

subplot(3,10,20);

for aa = 1:4
    plot(freq_v, data_2_plot_2{p,1}(2,:,aa),'linewidth',1,'color',[color_1])
    hold on
end

    plot(freq_v, data_2_plot_2{p,2}(2,:),'linewidth',3,'color',color_)

%set(gca, 'YTick', [])
xlim([3 12])
%%ylim([0 0.1])
title({'\fontsize{14}Averaged Trials',[]})



% dHPC
for ii = 1:t

    subplot(3,10,ii+20);

    for aa = 1:4
        plot(freq_v,data_2_plot_1{p,1}(3,:,ii,aa),'linewidth',1,'color',[color_1])    
        hold on
    end

    plot(freq_v,data_2_plot_1{p,2}(3,:,ii),'linewidth',3,'color',color_)


    %set(gca, 'YTick', [],'YColor','none')
    box off
    xlim([3 12])
    %ylim([-0.05 3.5])
    title({['\fontsize{14}CS-Trial Blk  ' num2str(ii)],[]})
end

subplot(3,10,30);

for aa = 1:4
    plot(freq_v, data_2_plot_2{p,1}(3,:,aa),'linewidth',1,'color',[color_1])
    hold on
end

    plot(freq_v, data_2_plot_2{p,2}(3,:),'linewidth',3,'color',color_)

%set(gca, 'YTick', [])
xlim([3 12])
%%ylim([0 0.1])
title({'\fontsize{14}Averaged Trials',[]})

%  h = findall(gcf, 'Type', 'axes');
% 
% for ii = 1:length(h)
%     ylim(h(ii), [0 1]);
% end


clear ('ii',"data_2_plot_2","data_2_plot_1")


%%


end

%% Save
% newStr1 = id(1:end-8);
% name_1  = strcat(Path,'/',newStr1,'_Granger_Spec_Mean_All_3');
% Path    = files.FilesLoaded{1,1}(ms).folder;
% saveas(gcf,name_1,'png') % save figure
%
% close all
%
% clear('ii','jj')
% clear('Combinations_2','Combinations_1','sub_idx','steps','data2plot_mean','cc','session','name_1','newStr1','path','data2plot')
