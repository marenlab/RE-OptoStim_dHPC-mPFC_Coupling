
% Filtering and getting angles from Hilbert transform

not = 7;
filters = [7.5 8.5];
fs = 1000;

% chR2
for ss = 1:size(all_data_chR2,1)

    for ii = 1:size(all_data_chR2{ss,1}{not,1},1)
        for tt = 1:size(all_data_chR2{ss,1}{not,1},3)
            for jj = 1:size(filters,1)

                temp = eegfilt(all_data_chR2{ss,1}{not,1}(ii,:,tt),fs,filters(jj,1),filters(jj,2));
                data_angle{ss,1}(ii,:,tt) = angle(hilbert(temp));

            end
        end
    end

end


% mcherry
for ss = 1:size(all_data_mCherry,1)

    for ii = 1:size(all_data_mCherry{ss,1}{not,1},1)
        for tt = 1:size(all_data_mCherry{ss,1}{not,1},3)
            for jj = 1:size(filters,1)

                temp = eegfilt(all_data_mCherry{ss,1}{not,1}(ii,:,tt),fs,filters(jj,1),filters(jj,2));
                data_angle{ss,2}(ii,:,tt) = angle(hilbert(temp));

            end
        end
    end

end

clear ('ss','ii','tt','jj','temp','not','filters')

%%

phase = reshape(data_angle, 3, 10000, 5, []);
phase = permute(data_angle, [1, 3, 2, 4]);
phase = reshape(data_angle, 3, 5*10000, 9);

phase = circ_mean(phase(:,:,1:5),[],3);
%%

nsurrog = 1000; % number of rearrangements
trials = 1:5; % select trials

% plv_null = zeros(2,nsurrog);
% phase_diff_shuffled_all = zeros(2,nsurrog);  % store permutations



for tt = 1:size(trials,2)

    for ss = 1:size(data_angle,1)
        for gg = 1:size(data_angle,2)
    
            % Baseline
            phase1 = reshape(data_angle{ss,gg}(1,:,trials(tt)), 1, []);
            phase2 = reshape(data_angle{ss,gg}(2,:,trials(tt)), 1, []);
            phase3 = reshape(data_angle{ss,gg}(3,:,trials(tt)), 1, []);

            % Baseline
            % phase1 = reshape(data_angle{ss,gg}(1,:,[trials(tt):trials(tt)+4]), 1, []);
            % phase2 = reshape(data_angle{ss,gg}(2,:,[trials(tt):trials(tt)+4]), 1, []);
            % phase3 = reshape(data_angle{ss,gg}(3,:,[trials(tt):trials(tt)+4]), 1, []);


            for ii = 1:nsurrog

                % shift_amount = randi([1, length(phase1)-1]);
                % shuffled_phase1 = circshift(phase1,shift_amount,2);  % permute signal phase1
                

                % dHPC - mPFC IL
                shuffled_phase1 = shuffle_esc(unwrap(phase1),fs,200)'; % modified function for Shift array circularly
                phase_diff_shuffled = angle(exp(1i * (unwrap(phase3) - shuffled_phase1)));      
                phase_diff_shuffled_all{ss,gg}(1,ii,tt) = circ_mean(phase_diff_shuffled,[],2);  

                z_shuff = mean(exp(1i * (unwrap(phase3) - shuffled_phase1)));  
                plv_null{ss,gg}(1,ii,tt) = abs(z_shuff)';                     

                % dHPC - mPFC PL
                shuffled_phase2 = shuffle_esc(unwrap(phase1),fs,200)'; % modified function for Shift array circularly
                phase_diff_shuffled = angle(exp(1i * (unwrap(phase3) - shuffled_phase2)));      
                phase_diff_shuffled_all{ss,gg}(2,ii,tt) = circ_mean(phase_diff_shuffled,[],2);  

                z_shuff = mean(exp(1i * (unwrap(phase3) - shuffled_phase2)));  
                plv_null{ss,gg}(2,ii,tt) = abs(z_shuff)';                      



            end
        end
    end
end





for ii = 1:size(plv_null,2)
        plv_z{5,ii} = cat(4,plv_z{:, ii});
        plv_z{5,ii} = median(plv_z{5,ii},4);

        plv_z{5,ii} = cat(4,plv_z{:, ii});
        plv_z{5,ii} = median(plv_z{5,ii},4);
        % phase_diff_shuffled_all{5,ii} = cat(4,phase_diff_shuffled_all{:,ii});
        % phase_diff_shuffled_all{5,ii} = median(phase_diff_shuffled_all{5,ii},4);

end





clear ('ii','tt','ss','gg','phase1','phase2','phase3','shuffled_phase1','shuffled_phase2','phase_diff_shuffled','z_shuff','nsurrog','trials')


%% Comparing z values

for ii = 1:size(plv_null,1)-1
    for jj = 1:size(plv_null,2)
        plv_z{ii,jj} = (plv_real{ii,jj} - squeeze(mean(plv_null{ii,jj},2)))./squeeze(std(plv_null{ii,jj},[],2));
    end
end

for ii = 1:size(plv_null,1)
    for jj = 1:size(plv_null,2)
        plv_p{ii,jj} = erfc(plv_z{ii,jj} / sqrt(2)) / 2;
    end
end


%% Plots

blk_trial = 9;

group = 2; % 1 chR2 / 2 mCherry

if group == 1
    figure
    set(gcf,'color','w');
    set(gcf, 'Position', [1522, -330, 1925, 680]);

    color_ = [.2, .6, 1]; % chR2
else
    color_ = [.8, 0, 0]; % mCherry
    hold on

end

hold on
titles_ = {'Animal 1','Animal 2','Animal 3','Animal 4'};

sgtitle({'\fontsize{18} \bf Extinction - 9th trials blk - 8 Hz',[]})

% dHPC <-> mPFC-IL
for cc = 1:size(plv_null,1)

    subplot(3,5,cc)
    h = histogram(normalize(plv_null{cc, group}(1,:,blk_trial),2),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)
    ylabel({'\fontsize{14} dHPC <-> mPFC-IL','\fontsize{12} pdf'})
    box off


    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end
    
    hold on
    plot([plv_z{cc, group}(1,blk_trial) plv_z{cc, group}(1,blk_trial)],[0 .1],'color',color_, 'linew',4);
    xlabel(['p = ' num2str(plv_p{cc, 1}(1,blk_trial))])
    ylabel({'\fontsize{14} mPFC-IL','\fontsize{12} pdf'})
    box off

    if cc~=5
        title({['\fontsize{14}' titles_{cc}],[]})
    else
        title({['\fontsize{14}' 'all animals'],[]})
    end

end

% dHPC <-> mPFC-PL
for cc = 1:size(plv_null,1)

    subplot(3,5,cc+10)
    h = histogram(normalize(plv_null{cc, group}(2,:,blk_trial),2),50,'Normalization','pdf');
    h.FaceColor = color_;
    h.FaceAlpha = .5;
    h.EdgeColor = 'none';
    xlim([-4 30])
    ylim([0 .5])
    xline(1.645)
    ylabel({'\fontsize{14} dHPC <-> mPFC-IL','\fontsize{12} pdf'})
    box off


    if cc~=1
        set(gca, 'YTick', [],'YColor','none')
    end
    
    hold on
    plot([plv_z{cc, group}(2,blk_trial) plv_z{cc, group}(2,blk_trial)],[0 .1],'color',color_, 'linew',4);
    xlabel(['p = ' num2str(plv_p{cc, 1}(2,blk_trial))])
    ylabel({'\fontsize{14} mPFC-IL','\fontsize{12} pdf'})
    box off

    if cc~=5
        title({['\fontsize{14}' titles_{cc}],[]})
    else
        title({['\fontsize{14}' 'all animals'],[]})
    end

end

%% Save

newStr = 'permut_chR2_mcherry_blk_';
path = '/Users/flavio/Desktop';
name = strcat(path,'/',newStr,num2str(blk_trial),'_PLV_extinction');

saveas(gcf,name,'png')
exportgraphics(gcf,strcat(name,'.eps'),'ContentType','vector')

close all

clear('cc','color_1','h','titles_','blk_trial')

%%

