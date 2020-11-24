if contains(d_v,'dorsal')
    if contains(to_use,'pre')
        temp_preD = temporal_sim
        rand_preD = temporal_unrelated
    elseif contains(to_use,'post')
        temp_postD = temporal_sim
        rand_postD = temporal_unrelated
    end
else contains(d_v,'ventral')
    if contains(to_use,'pre')
        temp_preV = temporal_sim
        rand_preV = temporal_unrelated
    elseif contains(to_use,'post')
        temp_postV = temporal_sim
        rand_postV = temporal_unrelated
    end
end

if exist('temp_preD') & exist('temp_postD') & exist('temp_preV') & exist('temp_postV')
    display('Proceed')
else
    display('Incomplete Variables')
end
%%
map1 = cbrewer('seq','Blues',2)
map2 = cbrewer('seq','Oranges',2)
map = [map1([1:2],:);map2]

figure('Position',[0 100 900 500]); hold on
% Dorsal Pre Unrelated
%bar(1,nanmean(temporal_unrelated),'FaceColor',[0 0 0])
scatter(ones(length(rand_preD),1)-.5,rand_preD,800,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0], 'MarkerFaceAlpha', .05 , 'MarkerEdgeAlpha', .15)
errorbar(.5,mean(rand_preD),std(rand_preD)/sqrt(numel(rand_preD)),'r','vertical','LineWidth',3)

% Dorsal Pre Temporal
%bar(2,nanmean(temporal_sim),'FaceColor',map(1,:))
scatter(ones(length(temp_preD),1),temp_preD,800,'MarkerFaceColor',map(1,:),'MarkerEdgeColor',map(1,:)*.8, 'MarkerFaceAlpha', .25 , 'MarkerEdgeAlpha', .75)
errorbar(1,mean(temp_preD),std(temp_preD)/sqrt(numel(temp_preD)),'r','vertical','LineWidth',3)


% Ventral Pre Unrelated
%bar(3,nanmean(hemi_sim),'FaceColor',map(2,:))
scatter(ones(length(rand_preV),1)+1.5,rand_preV,800,'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerFaceAlpha', .25 , 'MarkerEdgeAlpha', .75)
errorbar(2.5,mean(rand_preV),std(rand_preV)/sqrt(numel(rand_preV)),'r','vertical','LineWidth',3)

% Ventral Pre Temporal
%bar(4,nanmean(temp_cohort_sim),'FaceColor',map(3,:))
scatter(ones(length(temp_preV),1)+2,temp_preV,800,'MarkerFaceColor',map(2,:),'MarkerEdgeColor',map(2,:)*.8, 'MarkerFaceAlpha', .25 , 'MarkerEdgeAlpha', .75)
errorbar(3,nanmean(temp_preV),nanstd(temp_preV)/sqrt(numel(temp_preV)),'r','vertical','LineWidth',3)



% Dorsal Post Random
%bar(4,nanmean(temp_cohort_sim),'FaceColor',map(3,:))
scatter(ones(length(rand_postD),1)+3.5,rand_postD,800,'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerFaceAlpha', .25 , 'MarkerEdgeAlpha', .75)
errorbar(4.5,nanmean(rand_postD),nanstd(rand_postD)/sqrt(numel(rand_postD)),'r','vertical','LineWidth',3)

% Dorsal Post Temporal
%bar(4,nanmean(temp_cohort_sim),'FaceColor',map(3,:))
scatter(ones(length(temp_postD),1)+4,temp_postD,800,'MarkerFaceColor',map(3,:),'MarkerEdgeColor',map(3,:)*.8, 'MarkerFaceAlpha', .25 , 'MarkerEdgeAlpha', .75)
errorbar(5,nanmean(temp_postD),nanstd(temp_postD)/sqrt(numel(temp_postD)),'r','vertical','LineWidth',3)


% Ventral Post Random
%bar(4,nanmean(temp_cohort_sim),'FaceColor',map(3,:))
scatter(ones(length(rand_postV),1)+5.5,rand_postV,800,'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerFaceAlpha', .25 , 'MarkerEdgeAlpha', .75)
errorbar(6.5,nanmean(rand_postV),nanstd(rand_postV)/sqrt(numel(rand_postV)),'r','vertical','LineWidth',3)

% Ventral Post Temporal
%bar(4,nanmean(temp_cohort_sim),'FaceColor',map(3,:))
scatter(ones(length(temp_postV),1)+6,temp_postV,800,'MarkerFaceColor',map(4,:),'MarkerEdgeColor',map(4,:)*.8, 'MarkerFaceAlpha', .25 , 'MarkerEdgeAlpha', .75)
errorbar(7,nanmean(temp_postV),nanstd(temp_postV)/sqrt(numel(temp_postV)),'r','vertical','LineWidth',3)




xticks([.5 1 2.5 3 4.5 5 6.5 7])
xlim([0 7.5])
ylim([0 1])
xticklabels({'Random','Temporal','Random','Temporal','Random','Temporal','Random','Temporal'})
ylabel('Synapse Similarity')
set(gca,'FontSize',18)

pvals.postV = ranksum(temp_postV,rand_postV)
pvals.preV = ranksum(temp_preV,rand_preV)
pvals.postD = ranksum(temp_postD,rand_postD)
pvals.preD = ranksum(temp_preD,rand_preD)
