% PW permutation stim period

% all_data_chR2.mat and all_data_mCherr.mat need to be loaded before

% all_data_XXX
% First column: raw recordings
% Second column: 3–5 Hz calculation every 5 trials
% Third column: 8 Hz calculation every 5 trials

% null_dist
% Within each group:
% First column: null distribution
% Second column: z-value
% Third column: p-value (one-tailed)

% Careful: note that in the original variable, the 3–5 Hz data is in the third column, and the 8 Hz data is in the second column.

%%

window = 2000;
noverlap = 1000;
nfft = 2^15;
fs = 1000;

n_iter = 1000;

low_freq  = 3;
high_freq = 5;

not = 7;

scale = 1/0.153; % to fiz plexon scale

numshf  = 1000; % number of shuffled segments
nsurrog = 1000; % number of rearrangements

trials = 1:5:45;


% mCherry

for tt = 1:size(trials,2)


    for ss = 1:size(all_data_mCherry,1)

        full_signal = reshape(all_data_mCherry{ss,1}{not,1}(:,:,[trials(tt):trials(tt)+4]), 3, []);

        for ch = 1:3
            for ii = 1:nsurrog
                shuffled(ch,:) = shuffle_esc(full_signal(ch,:),fs,numshf);
                [pxx(ch,:), f] = pwelch(shuffled(ch,:), window, noverlap, nfft, fs);

                idx_band_norm = (f >= 2) & (f <= 12);
                pxx(ch,:) = scale.*(pxx(ch,:)./sum(pxx(ch,idx_band_norm)));

                idx_band_1 = (f >= 3) & (f <= 5);
                null_dist.mCherry.f_3_5Hz{ss,1}(ch,ii,tt) = mean(pxx(ch,idx_band_1));

                idx_band_2 = (f >= 7.8) & (f <= 8.2);
                null_dist.mCherry.f_7_9Hz{ss,1}(ch,ii,tt) = mean(pxx(ch,idx_band_2));
            end
        end

        null_dist.mCherry.f_3_5Hz_Z{ss,1}(:,:,tt) = normalize(null_dist.mCherry.f_3_5Hz{ss,1}(:,:,tt),2);        
        null_dist.mCherry.f_7_9Hz_Z{ss,1}(:,:,tt) = normalize(null_dist.mCherry.f_7_9Hz{ss,1}(:,:,tt),2);

    end

end


% chR2

for tt = 1:size(trials,2)


    for ss = 1:size(all_data_chR2,1)

        full_signal = reshape(all_data_chR2{ss,1}{not,1}(:,:,[trials(tt):trials(tt)+4]), 3, []);

        for ch = 1:3
            for ii = 1:nsurrog
                shuffled(ch,:) = shuffle_esc(full_signal(ch,:),fs,numshf);
                [pxx(ch,:), f] = pwelch(shuffled(ch,:), window, noverlap, nfft, fs);

                idx_band_norm = (f >= 2) & (f <= 12);
                pxx(ch,:) = scale.*(pxx(ch,:)./sum(pxx(ch,idx_band_norm)));

                idx_band_1 = (f >= 3) & (f <= 5);
                null_dist.chR2.f_3_5Hz{ss,1}(ch,ii,tt) = mean(pxx(ch,idx_band_1));

                idx_band_2 = (f >= 7.8) & (f <= 8.2);
                null_dist.chR2.f_7_9Hz{ss,1}(ch,ii,tt) = mean(pxx(ch,idx_band_2));
            end
        end

        null_dist.chR2.f_3_5Hz_Z{ss,1}(:,:,tt) = normalize(null_dist.chR2.f_3_5Hz{ss,1}(:,:,tt),2);
        null_dist.chR2.f_7_9Hz_Z{ss,1}(:,:,tt) = normalize(null_dist.chR2.f_7_9Hz{ss,1}(:,:,tt),2);

    end

end

clear('idx_band_2', 'nn','trials','tt', 'full_signal', 'n_iter', 'window', 'noverlap', 'nfft', 'fs', 'ii','idx_band_1', 'f', 'low_freq', 'high_freq', 'ss', 'numshf', 'nsurrog', 'not', 'scale', 'ch', 'shuffled','pxx','idx_band_norm')

%% Comparing z values


for nn = 1:size(null_dist.mCherry.f_3_5Hz_Z,1)

    null_dist.mCherry.f_3_5Hz_Z{nn,2} = (all_data_mCherry{nn,3} - squeeze(mean(null_dist.mCherry.f_3_5Hz{nn,1},2))) ./ squeeze(std(null_dist.mCherry.f_3_5Hz{nn,1},[],2));   
    null_dist.mCherry.f_3_5Hz_Z{nn,3} = erfc(null_dist.mCherry.f_3_5Hz_Z{nn,2} / sqrt(2)) / 2;

end


for nn = 1:size(null_dist.mCherry.f_7_9Hz_Z,1)

    null_dist.mCherry.f_7_9Hz_Z{nn,2} = (all_data_mCherry{nn,2} - squeeze(mean(null_dist.mCherry.f_7_9Hz{nn,1},2))) ./ squeeze(std(null_dist.mCherry.f_7_9Hz{nn,1},[],2));   
    null_dist.mCherry.f_7_9Hz_Z{nn,3} = erfc(null_dist.mCherry.f_7_9Hz_Z{nn,2} / sqrt(2)) / 2;

end





for nn = 1:size(null_dist.chR2.f_3_5Hz_Z,1)

    null_dist.chR2.f_3_5Hz_Z{nn,2} = (all_data_chR2{nn,3} - squeeze(mean(null_dist.chR2.f_3_5Hz{nn,1},2))) ./ squeeze(std(null_dist.chR2.f_3_5Hz{nn,1},[],2));   
    null_dist.chR2.f_3_5Hz_Z{nn,3} = erfc(null_dist.chR2.f_3_5Hz_Z{nn,2} / sqrt(2)) / 2;

end


for nn = 1:size(null_dist.chR2.f_7_9Hz_Z,1)

    null_dist.chR2.f_7_9Hz_Z{nn,2} = (all_data_chR2{nn,2} - squeeze(mean(null_dist.chR2.f_7_9Hz{nn,1},2))) ./ squeeze(std(null_dist.chR2.f_7_9Hz{nn,1},[],2));   
    null_dist.chR2.f_7_9Hz_Z{nn,3} = erfc(null_dist.chR2.f_7_9Hz_Z{nn,2} / sqrt(2)) / 2;

end




clear('nn','freq_to_')

%% Plots

data_to_plot_ = null_dist;
color_ = [.2, .6, 1]; % chR2

blk_trial = 1;

% data_to_plot_ = null_dist.chR2.f_3_5Hz_Z;
% color_ = [.2, .6, 1]; % chR2

% data_to_plot_ = null_dist.chR2.f_7_9Hz_Z;
% color_ = [.2, .6, 1]; % chR2




for ii = 1:size(data_to_plot_,2)
    if ii == 1
        data_to_plot_{5,ii} = cat(4,data_to_plot_{:, ii});
        data_to_plot_{5,ii} = reshape(permute(data_to_plot_{5,ii}, [1 2 4 3]), 3, 1000*4, 9);
    else
        data_to_plot_{5,ii} = median(cat(3,data_to_plot_{:, ii}),3);
    end
end



figure
set(gcf,'color','w');
set(gcf, 'Position', [1522, -330, 1925, 680]);


hold on
titles_ = {'Animal 1','Animal 2','Animal 3','Animal 4'};

sgtitle({'\fontsize{18} \bf Extinction - 2th trials blk - 8 Hz',[]})

%IL
for cc = 1:size(data_to_plot_,1)

    subplot(3,5,cc)
    h = histogram(data_to_plot_{cc, 1}(1,:,blk_trial),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)

    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end
    
    hold on
    plot([data_to_plot_{cc, 2}(1,blk_trial) data_to_plot_{cc, 2}(1,blk_trial)],[0 .1],'color',color_, 'linew',4);
    xlabel(['p = ' num2str(data_to_plot_{cc, 3}(1,blk_trial))])
    ylabel({'\fontsize{14} mPFC-IL','\fontsize{12} pdf'})
    box off

    if cc~=5
        title({['\fontsize{14}' titles_{cc}],[]})
    else
        title({['\fontsize{14}' 'all animals'],[]})
    end

end



%PL
for cc = 1:size(data_to_plot_,1)

    subplot(3,5,cc+5)
    h = histogram(data_to_plot_{cc, 1}(2,:,blk_trial),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)

    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end

    hold on
    plot([data_to_plot_{cc, 2}(2,blk_trial) data_to_plot_{cc, 2}(2,blk_trial)],[0 .1],'color',color_,'linew',4);
    xlabel(['p = ' num2str(data_to_plot_{cc, 3}(2,blk_trial))])
    ylabel({'\fontsize{14} mPFC-PL','\fontsize{12} pdf'})
    box off


end


%dHPC
for cc = 1:size(data_to_plot_,1)

    subplot(3,5,cc+10)
    h = histogram(data_to_plot_{cc, 1}(3,:,blk_trial),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)

    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end

    hold on
    plot([data_to_plot_{cc, 2}(3,blk_trial) data_to_plot_{cc, 2}(3,blk_trial)],[0 .1],'color',color_,'linew',4);
    xlabel(['p = ' num2str(data_to_plot_{cc, 3}(3,blk_trial))])
    ylabel({'\fontsize{14} dHPC','\fontsize{12} pdf'})
    box off


end

legend('null distribuition','observed value')
legend boxoff  
clear('cc','color_','h','titles_')






% data_to_plot_ = null_dist.mCherry.f_3_5Hz_Z;
% color_ = [.8, 0, 0]; % mCherry

data_to_plot_ = null_dist.mCherry.f_7_9Hz_Z;
color_ = [.8, 0, 0]; % mCherry

for ii = 1:size(data_to_plot_,2)
    if ii == 1
        data_to_plot_{5,ii} = cat(4,data_to_plot_{:, ii});
        data_to_plot_{5,ii} = reshape(permute(data_to_plot_{5,ii}, [1 2 4 3]), 3, 1000*4, 9);
    else
        data_to_plot_{5,ii} = mean(cat(3,data_to_plot_{:, ii}),3);
    end
end



titles_ = {'Animal 1','Animal 2','Animal 3','Animal 4'};
%IL
for cc = 1:size(data_to_plot_,1)

    subplot(3,5,cc)
    h = histogram(data_to_plot_{cc, 1}(1,:,blk_trial),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)

    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end
    
    hold on
    plot([data_to_plot_{cc, 2}(1,blk_trial) data_to_plot_{cc, 2}(1,blk_trial)],[0 .1],'color',color_, 'linew',4);
    xlabel(['p = ' num2str(data_to_plot_{cc, 3}(1,blk_trial))])
    ylabel({'\fontsize{14} mPFC-IL','\fontsize{12} pdf'})
    box off

    if cc~=5
        title({['\fontsize{14}' titles_{cc}],[]})
    else
        title({['\fontsize{14}' 'all animals'],[]})
    end

end



%PL
for cc = 1:size(data_to_plot_,1)

    subplot(3,5,cc+5)
    h = histogram(data_to_plot_{cc, 1}(2,:,blk_trial),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)

    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end

    hold on
    plot([data_to_plot_{cc, 2}(2,blk_trial) data_to_plot_{cc, 2}(2,blk_trial)],[0 .1],'color',color_,'linew',4);
    xlabel(['p = ' num2str(data_to_plot_{cc, 3}(2,blk_trial))])
    ylabel({'\fontsize{14} mPFC-PL','\fontsize{12} pdf'})
    box off


end


%dHPC
for cc = 1:size(data_to_plot_,1)

    subplot(3,5,cc+10)
    h = histogram(data_to_plot_{cc, 1}(3,:,blk_trial),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)

    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end

    hold on
    plot([data_to_plot_{cc, 2}(3,blk_trial) data_to_plot_{cc, 2}(3,blk_trial)],[0 .1],'color',color_,'linew',4);
    xlabel(['p = ' num2str(data_to_plot_{cc, 3}(3,blk_trial))])
    ylabel({'\fontsize{14} dHPC','\fontsize{12} pdf'})
    box off


end

legend('null distribuition','observed value')
legend boxoff  
clear('cc','color_','h')
