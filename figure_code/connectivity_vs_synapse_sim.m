to_test = 'output'


if contains(to_test,'output')
    map = cbrewer('seq','Blues',4)
elseif contains(to_test,'input')
    map = cbrewer('seq','Oranges',4)

end


%% Connectivity Similarity
figure; hold on
Violin(connectivity_data.Random,1,'ViolinColor',[0,0,0],'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])
Violin(connectivity_data.Temporal,2,'ViolinColor',map(1,:),'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])
Violin(connectivity_data.Hemilineage,3,'ViolinColor',map(2,:),'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])
Violin(connectivity_data.Temporal_Cohort,4,'ViolinColor',map(3,:),'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])
xticks([1 2 3 4])
xlim([0 5])
ylim([0 1])
xticklabels({'Random','Temporal','Hemilineage','HL-Temporal Cohort'})
ylabel('Input Connectivity Similarity')
set(gca,'FontSize',18)
%% Stats
data = 0 
group = 0
data = connectivity_data.Random;
group = ones(numel(connectivity_data.Random),1);
data = [data; connectivity_data.Temporal]
group = [group ; ones(numel(connectivity_data.Temporal),1) + 1];
data = [data; connectivity_data.Hemilineage];
group = [group ; ones(numel(connectivity_data.Hemilineage),1) + 2];
data = [data ; connectivity_data.Temporal_Cohort];
group = [group ; ones(numel(connectivity_data.Temporal_Cohort),1) + 3];


% Anova / multicple comparisons
group(isnan(data)) = []
data(isnan(data)) = []
[p,tbl,stats] = anovan(data,{group})
[results,means] = multcompare(stats,'CType','bonferroni')
%%
to_remove = isnan(connectivity_data.Random);
connectivity_data.Random(to_remove) = 0; clear to_remove

to_remove = isnan(connectivity_data.Temporal);
connectivity_data.Temporal(to_remove) = 0; clear to_remove

to_remove = isnan(connectivity_data.Hemilineage);
connectivity_data.Hemilineage(to_remove) = 0; clear to_remove

to_remove = isnan(connectivity_data.Temporal_Cohort);
connectivity_data.Temporal_Cohort(to_remove) = 0; clear to_remove


figure; subplot(1,4,1)
%scatter(similarity_data.Random,connectivity_data.Random,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
lm_rand = fitlm(similarity_data.Random,connectivity_data.Random)
h1 = plot(lm_rand)
h1(1).Color = 'k'
h1(1).Marker = '.'
h1(1).MarkerSize = 30
h1(2).Color = 'k'
h1(3).Color = 'k'
h1(4).Color = 'k'
xlim([0 1])
ylim([0 1])

subplot(1,4,2)
%scatter(similarity_data.Temporal,connectivity_data.Temporal,'MarkerFaceColor',map(1,:),'MarkerEdgeColor',map(1,:),'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
lm_temp = fitlm(similarity_data.Temporal,connectivity_data.Temporal)
h2 = plot(lm_temp)
h2(1).Color = map(1,:)
h2(1).Marker = '.'
h2(1).MarkerSize = 30
h2(2).Color = map(1,:)
h2(3).Color = map(1,:)
h2(4).Color = map(1,:)
xlim([0 1])
ylim([0 1])

subplot(1,4,3)
%scatter(similarity_data.Hemilineage,connectivity_data.Hemilineage,'MarkerFaceColor',map(2,:),'MarkerEdgeColor',map(2,:),'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
lm_hemi = fitlm(similarity_data.Hemilineage,connectivity_data.Hemilineage)
h3 = plot(lm_hemi)
h3(1).Color = map(2,:)
h3(1).Marker = '.'
h3(1).MarkerSize = 30
h3(2).Color = map(2,:)
h3(3).Color = map(2,:)
h3(4).Color = map(2,:)
xlim([0 1])
ylim([0 1])

subplot(1,4,4)
%scatter(similarity_data.Temporal_Cohort,connectivity_data.Temporal_Cohort,'MarkerFaceColor',map(3,:),'MarkerEdgeColor',map(3,:),'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
lm_tc = fitlm(similarity_data.Temporal_Cohort,connectivity_data.Temporal_Cohort)
h4 = plot(lm_tc)
h4(1).Color = map(3,:)
h4(1).Marker = '.'
h4(1).MarkerSize = 30
h4(2).Color = map(3,:)
h4(3).Color = map(3,:)
h4(4).Color = map(3,:)
xlim([0 1])
ylim([0 1])
%%
sim = {similarity_data.Random,similarity_data.Temporal,similarity_data.Hemilineage,similarity_data.Temporal_Cohort}
con = {connectivity_data.Random,connectivity_data.Temporal,connectivity_data.Hemilineage,connectivity_data.Temporal_Cohort}
x = [.025:.05:1]

figure; hold on
clear mean_group
for i = 1:length(sim)
    y = sim{i}(isnan(sim{i})==0)
    [N,edges,bin] = histcounts(y,[0:.05:1])
     
     for j = 1:length(edges)-1
         group = con{i}(bin == j)
         mean_group(j) = mean(group)
%         %bar(x(j),nanmean(group),.05,'EdgeColor',map(i,:),'FaceAlpha',.001)
%         %errorbar(x(j),nanmean(group)/sqrt(numel(group)),nanstd(group),'k')
     end
     plot(x,mean_group,'Color',map(i,:),'LineWidth',3)
     

end
%%
bins = [.05:.1:1]
figure; subplot(3,1,1)
bar(bins,group_count)
ylabel('Frequency')


