%% Load dataset
load mark2020_neurons_temporal_cohorts.mat
%% Get lineage groupings and calculate synapse similarities

lineage_grouping = an_in.Lineage_Index(an_in.Side_Index == 0)

similarity.pre.left = synapse_similarity_v2(nl(an_in.Side_Index == 0),2000,3,[],1)
similarity.pre.right = synapse_similarity_v2(nl(an_in.Side_Index == 1),2000,3,[],1)
similarity.pre.combined = (similarity.pre.left + similarity.pre.right)*.5

similarity.post.left = synapse_similarity_v2(nl(an_in.Side_Index == 0),2000,3,[],2)
similarity.post.right = synapse_similarity_v2(nl(an_in.Side_Index == 1),2000,3,[],2)
similarity.post.combined = (similarity.post.left + similarity.post.right)*.5
%% presynapses

%For each lineage, get the similarity scores of intra lineage pairs vs
%inter lineage pairs.

for i = 1:max(lineage_grouping)
    related_index = lineage_grouping == i;
    unrelated_index = boolean((related_index*-1)+1);
    
    sim_mat_related = similarity.pre.combined(related_index,related_index);
    sim_mat_unrelated = similarity.pre.combined(related_index,unrelated_index);
    
    related_vals = sim_mat_related(boolean(triu(ones(size(sim_mat_related)),1)));
    unrelated_vals = sim_mat_unrelated(:);
    
    related{i} = related_vals;
    unrelated{i} = unrelated_vals; clear sim_mat_related and sim_mat_unrelated
end


intra_lineage_pre = cat(1,related{:});
inter_lineage_pre = cat(1,unrelated{:}); clear related and unrelated

real_mean_dif = mean(intra_lineage_pre) - mean(inter_lineage_pre)


%% postsynapses

%For each lineage, get the similarity scores of intra lineage pairs vs
%inter lineage pairs.
for i = 1:max(lineage_grouping)
    related_index = lineage_grouping == i;
    unrelated_index = boolean((related_index*-1)+1);
    
    sim_mat_related = similarity.post.combined(related_index,related_index);
    sim_mat_unrelated = similarity.post.combined(related_index,unrelated_index);
    
    related_vals = sim_mat_related(boolean(triu(ones(size(sim_mat_related)),1)));
    unrelated_vals = sim_mat_unrelated(:);
    
    related{i} = related_vals;
    unrelated{i} = unrelated_vals ; clear sim_mat_related and sim_mat_unrelated
end


intra_lineage_post = cat(1,related{:});
inter_lineage_post = cat(1,unrelated{:}); clear related and unrelated

real_mean_dif = mean(intra_lineage_post) - mean(inter_lineage_post)
%% Plot mean similarities


figure('pos',[0,0,400,800]); hold on
%scatter(ones(length(inter_lineage_pre),1),inter_lineage_pre,800,'MarkerFaceColor','k','MarkerFaceAlpha',.05,'MarkerEdgeColor','k','MarkerEdgeAlpha',.25)
%bar(1,mean(inter_lineage_pre),'FaceColor','k')
errorbar(1,mean(inter_lineage_pre),std(inter_lineage_pre)/sqrt(numel(inter_lineage_pre)),'k','LineWidth',2)
Violin(inter_lineage_pre,1,'ViolinColor',[0 0 0],'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])


%scatter(ones(length(intra_lineage_pre),1)+1,intra_lineage_pre,800,'MarkerFaceColor',[0.6196 0.7922 0.8824],'MarkerFaceAlpha',.15,'MarkerEdgeColor',[0.1922 0.5098 0.7412],'MarkerEdgeAlpha',.25)
%bar(2,mean(intra_lineage_pre),'FaceColor',[0.6196 0.7922 0.8824])
%errorbar(2,mean(intra_lineage_pre),std(intra_lineage_pre)/sqrt(numel(intra_lineage_pre)),'k','LineWidth',2)
Violin(intra_lineage_pre,2,'ViolinColor',[0.6196 0.7922 0.8824],'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])


%scatter(ones(length(intra_lineage_post),1)+3,intra_lineage_post,800,'MarkerFaceColor','k','MarkerFaceAlpha',.05,'MarkerEdgeColor','k','MarkerEdgeAlpha',.25)
%bar(4,mean(inter_lineage_post),'FaceColor','k')
%errorbar(4,mean(inter_lineage_post),std(inter_lineage_post)/sqrt(numel(inter_lineage_post)),'k','LineWidth',2)
Violin(inter_lineage_post,4,'ViolinColor',[0 0 0],'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])


%scatter(ones(length(inter_lineage_post),1)+4,inter_lineage_post,800,'MarkerFaceColor',[0.9922 0.6824 0.4196],'MarkerFaceAlpha',.05,'MarkerEdgeColor',[0.9020 0.3333 0.0510],'MarkerEdgeAlpha',.25)
%bar(5,mean(intra_lineage_post),'FaceColor',[0.9922 0.6824 0.4196])
%errorbar(5,mean(intra_lineage_post),std(intra_lineage_post)/sqrt(numel(intra_lineage_post)),'k','LineWidth',2)
Violin(intra_lineage_post,5,'ViolinColor',[0.9922 0.6824 0.4196],'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])



%ylim([0 1])
xticks([1 2 4 5])
xticklabels({'Inter Lineage','Intra Lineage','Inter Lineage','Intra Lineage'})
xlim([.5 5.5])
set(gca,'FontSize',18)
ylabel('Synapse Similarity')

[p_pre] = ranksum(inter_lineage_pre,intra_lineage_pre)
[p_post] = ranksum(inter_lineage_post,intra_lineage_post)

%% Empirical distributions 

[f_pre_intra,x_pre_intra] = ecdf(intra_lineage_pre)
[f_pre_inter,x_pre_inter] =  ecdf(inter_lineage_pre)

[f_post_intra,x_post_intra] = ecdf(intra_lineage_post)
[f_post_inter,x_post_inter] =  ecdf(inter_lineage_post)


figure; hold on;
plot(x_pre_intra,f_pre_intra,'Color',[0.6196 0.7922 0.8824],'LineWidth',3)
plot(x_pre_inter,f_pre_inter,'Color',[0.6196 0.7922 0.8824]*.95,'LineWidth',3,'LineStyle','--')

plot(x_post_intra,f_post_intra,'Color',[0.9922 0.6824 0.4196],'LineWidth',3)
plot(x_post_inter,f_post_inter,'Color',[0.9922 0.6824 0.4196]*.95,'LineWidth',3,'LineStyle','--')

xlabel('Synapse Similarity')
ylabel('Frequency')
set(gca,'FontSize',24)

%% Permutation test
for j = 1:10000
    
    randomized_grouping = lineage_grouping(randperm(length(lineage_grouping)));
    
    for i = 1:max(randomized_grouping)
        related_index = randomized_grouping == i;
        unrelated_index = boolean((related_index*-1)+1);

        sim_mat_related = similarity.combined(related_index,related_index);
        sim_mat_unrelated = similarity.combined(related_index,unrelated_index);

        related_vals = sim_mat_related(boolean(triu(ones(size(sim_mat_related)),1)));
        unrelated_vals = sim_mat_unrelated(:);

        related{i} = related_vals;
        unrelated{i} = unrelated_vals;
    end
    
    intra_random = cat(1,related{:});
    inter_random = cat(1,unrelated{:}); clear related and unrelated
    
    random_mean_dif(j) = median(intra_random) - median(inter_random);
    
end
    
figure; hold on
histogram(random_mean_dif)

%% ANOSIM analysis 

stats.lineage = f_anosim(1-similarity.combined,lineage_grouping,[],10000,0,1)

an_half = an_in(1:2:end,:)
unique_temporal = unique(an_half(:,[2,4,6]),'rows')
[~,temporal_grouping] = ismember(an_half(:,[2,4,6]),unique_temporal,'rows')
stats.temporal = f_anosim(1-similarity.combined,temporal_grouping,[],10000,0,1)
%%
unique_hemilineage = unique(an_half(:,[2,4]),'rows')
[~,hemilineage_grouping] = ismember(an_half(:,[2,4]),unique_hemilineage,'rows')
stats.hemilineage = f_anosim(1-similarity.combined,hemilineage_grouping,[],10000,0,1)

%% dv grouping
[~,dv_grouping] = ismember(an_half.DV_Index,[0,1])
stats.dv = f_anosim(1-similarity.combined,dv_grouping,[],10000,0,1)

%% P-vals 
clear temporal_cohort_index and hl_mat

for i = 1:max(hemilineage_grouping)
    hemilineage = i

    temporal_cohort_index = temporal_grouping(hemilineage_grouping == hemilineage)
    hl_mat = similarity.combined(hemilineage_grouping == hemilineage,hemilineage_grouping == hemilineage)

    result.p = f_anosim(1-hl_mat,temporal_cohort_index,0,10000,0,1)
    p(i) = result.p
end