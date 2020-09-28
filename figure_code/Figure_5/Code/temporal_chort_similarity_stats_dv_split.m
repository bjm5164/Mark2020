load mark2020_neurons_temporal_cohorts.mat

% Remove NB1-2D from this analysis.
NB12D_index = boolean(strcmp(an_in.Lineage,'NB1-2') & an_in.DV_Index(:) == 1);
nl(NB12D_index) = [];
an_in(NB12D_index,:) = []; clear NB12D_index

%% Separate into dorsal and ventral hemilineages
dorsal_neurons = nl(an_in.DV_Index == 1)
an_dorsal = an_in(an_in.DV_Index == 1,:)
ventral_neurons = nl(an_in.DV_Index == 0)
an_ventral = an_in(an_in.DV_Index == 0,:)

%% Calculate similarity for either dorsal or ventral hemilineages
d_v = 'ventral'

if contains(d_v,'dorsal')
    nl_to_test = dorsal_neurons
    an_to_test = an_dorsal
elseif contains(d_v,'ventral')
    nl_to_test = ventral_neurons
    an_to_test = an_ventral
else errot('wrong d/v input')
end


presim_mat_l = synapse_similarity_v2(nl_to_test(an_to_test.Side_Index == 0),2000,3,[],1);
presim_mat_r = synapse_similarity_v2(nl_to_test(an_to_test.Side_Index == 1),2000,3,[],1);
presim = .5*(presim_mat_l+presim_mat_r)

postsim_mat_l = synapse_similarity_v2(nl_to_test(an_to_test.Side_Index == 0),2000,3,[],2);
postsim_mat_r = synapse_similarity_v2(nl_to_test(an_to_test.Side_Index == 1),2000,3,[],2);
postsim = .5*(postsim_mat_l+postsim_mat_r)

sim_mat = (presim+postsim)*.5
d_edge_total = [nl_to_test(an_to_test.Side_Index == 0).Temporal_Cohort];
an_half = an_to_test(1:2:end,:);
nl_half_skids = [nl_to_test(1:2:end).SkIDs];


figure; hold on; histogram(presim,[0:.1:1],'Normalization','probability'); histogram(postsim,[0:.1:1],'Normalization','probability')
%%

% Define which similarity matrix to use
%similarity_matrix = corrcov(presim_mat_l-presim_mat_r)
to_use = 'pre'

if contains(to_use,'pre')
    similarity_matrix = presim;
elseif contains(to_use,'post')
    similarity_matrix = postsim;
end

% Unrelated neurons: Neurons not related by lineage or birth order

unrelated_matrix = zeros(size(similarity_matrix))-1;
an_array = table2array(an_half(:,[2,4,6]));

clear nt_sim and t_sim and h_sim and nh_sim and tc_sim and ntc_sim
for i = 1:size(an_array,1)
    for j = 1:size(an_array,1)
        if numel(intersect(an_array(i,:),an_array(j,:))) > 1;
            unrelated_matrix(i,j) = NaN;
        else
            unrelated_matrix(i,j) = similarity_matrix(i,j);
        end
    end
end
sur = unrelated_matrix(boolean(triu(ones(size(unrelated_matrix)),1))) % Get pairwise comparisons of neurons in temporal cohort i and all neurons not in temporal cohort i
sim_unrelated = sur(:)
nt_sim{i} = sim_unrelated;


% Temporal
temporal_index = unique(an_half(:,6)); % Possible temporal cohorts

for i = 1:length(temporal_index.Temporal_Cohort)
    % Find neurons in temporal cohort i
    related_index = find(ismember(an_half(:,6),temporal_index(i,:))) 
    %get the pairwise similarities between them
    sr = similarity_matrix(related_index,related_index) 
    
    % Remove comparisons for temporal cohorts within hemilineages. 
    for ii = 1:length(related_index)
        for j = 1:length(related_index)
            if an_half.Lineage_Index(related_index(ii)) == an_half.Lineage_Index(related_index(j))
                sr(ii,j) = NaN
            end
        end
    end
    % Get only 1 side of the matrix to avoid self-self comparisons and
    % duplicates.    
    sr = sr(boolean(triu(ones(size(sr)),1)));
    sim_related = sr(:)
    t_sim{i} = sim_related; 
    clear sim_related and sim_unrelated and sr and sur
end

% Hemilineage 
hl_index = unique(an_half(:,[2,4]),'rows') % Possible hemilineages

for i = 1:length(hl_index.DV_Index)
    % Find neurons in hemilineage i
    related_index_h = find(ismember(an_half(:,[2,4]),hl_index(i,:),'rows')); 
    % get pairwise similarities between them
    sr_h = similarity_matrix(related_index_h,related_index_h); 
 
    % Remove comparisons for temporal cohorts within hemilineages. 
    for ii = 1:length(related_index_h)
        for j = 1:length(related_index_h)
            if an_half.Temporal_Cohort(related_index_h(ii)) == an_half.Temporal_Cohort(related_index_h(j))
                sr_h(ii,j) = NaN
            end
        end
    end
    
    % Get only 1 side of the matrix to avoid self-self comparisons and
    % duplicates. 
    sr_h = sr_h(boolean(triu(ones(size(sr_h)),1)))
    
    % Get pairwise comparisons of neurons in hemilineage i and all neurons not in hemilineage i
    sur_h = similarity_matrix(related_index_h,boolean((related_index_h*-1)+1)); 
    sim_related = sr_h(:);
    sim_unrelated = sur_h(:);
    h_sim{i} = sim_related; 
    nh_sim{i} = sim_unrelated; clear sim_related and sim_unrelated and sr_h and sur_h
end

% Temporal Cohort
tc_index = unique(an_half(:,[2,4,6]),'rows') % Possible temporal cohorts

for i = 1:length(tc_index.DV_Index)
    related_index = ismember(an_half(:,[2,4,6]),tc_index(i,:),'rows') % Find neurons in temporal cohort i
    sr = similarity_matrix(related_index,related_index)
    sr = sr(boolean(triu(ones(size(sr)),1)))
    sur = similarity_matrix(related_index,boolean((related_index*-1)+1));
    sim_unrelated = sur(:);
    sim_related = sr(:);
    %sim_related(sim_related == 1) = []
    tc_sim{i} = sim_related;
    ntc_sim{i} = sim_unrelated; clear sim_related and sim_unrelated and sr and sur
    

end


    

% 
% for i = 1:length(tc_sim)
%     if isempty(tc_sim{i})
%         idx(i) = 1
%     else
%         idx(i) = 0
%     end
% end
% 
% tc_sim(boolean(idx)) = []
%% Quantify similarity of groups

if contains(to_use,'pre')
    map = cbrewer('seq','Blues',3)
elseif contains(to_use,'post')
    map = cbrewer('seq','Oranges',3)

end


temporal_sim = cat(1,t_sim{:});
temporal_unrelated = cat(1,nt_sim{:});

hemi_sim = cat(1,h_sim{:});
hemi_unrelated = cat(1,nh_sim{:});

temp_cohort_sim = cat(1,tc_sim{:});
temp_cohort_unrelated = cat(1,ntc_sim{:});

temporal_unrelated(isnan(temporal_unrelated)) = []
temporal_sim(isnan(temporal_sim)) = []
hemi_sim(isnan(hemi_sim)) = []
temp_cohort_sim(isnan(temp_cohort_sim)) = []



% Stats
[p_temporal_vs_nontemporal,h,stat] = ranksum(temporal_sim,temporal_unrelated)
[p_hemilineage_vs_temporal,h,stat] = ranksum(temporal_sim,hemi_sim)
[p_hemilineage_vs_temporal_cohort,h,stat] = ranksum(hemi_sim,temp_cohort_sim)

figure('Position',[0 100 900 500]); hold on
% Unrelated neurons
%bar(1,nanmean(temporal_unrelated),'FaceColor',[0 0 0])
scatter(ones(length(temporal_unrelated),1),temporal_unrelated,800,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0], 'MarkerFaceAlpha', .05 , 'MarkerEdgeAlpha', .15)
errorbar(1,nanmean(temporal_unrelated),nanstd(temporal_unrelated)/sqrt(numel(temporal_unrelated)),'k','vertical','LineWidth',3)

% Temporal cohorts
%bar(2,nanmean(temporal_sim),'FaceColor',map(1,:))
scatter(ones(length(temporal_sim),1)+1,temporal_sim,800,'MarkerFaceColor',map(1,:),'MarkerEdgeColor',map(1,:)*.8, 'MarkerFaceAlpha', .25 , 'MarkerEdgeAlpha', .75)
errorbar(2,nanmean(temporal_sim),nanstd(temporal_sim)/sqrt(numel(temporal_sim)),'k','vertical','LineWidth',3)

%Hemilineages
%bar(3,nanmean(hemi_sim),'FaceColor',map(2,:))
scatter(ones(length(hemi_sim),1)+2,hemi_sim,800,'MarkerFaceColor',map(2,:),'MarkerEdgeColor',map(2,:)*.8, 'MarkerFaceAlpha', .25 , 'MarkerEdgeAlpha', .75)
errorbar(3,nanmean(hemi_sim),nanstd(hemi_sim)/sqrt(numel(hemi_sim)),'k','vertical','LineWidth',3)

% Temporal cohorts of hemilineages
%bar(4,nanmean(temp_cohort_sim),'FaceColor',map(3,:))
scatter(ones(length(temp_cohort_sim),1)+3,temp_cohort_sim,800,'MarkerFaceColor',map(3,:),'MarkerEdgeColor',map(3,:)*.8, 'MarkerFaceAlpha', .25 , 'MarkerEdgeAlpha', .75)
errorbar(4,nanmean(temp_cohort_sim),nanstd(temp_cohort_sim)/sqrt(numel(temp_cohort_sim)),'k','vertical','LineWidth',3)


xticks([1 2 3 4])
xlim([0 5])
ylim([0 1])
xticklabels({'Random','Temporal','Hemilineage','HL-Temporal Cohort'})
ylabel('Synapse Similarity')
set(gca,'FontSize',18)



%Cumulative histogram of similarity distributions. 
figure; hold on
tc_dist = histcounts(temp_cohort_sim,[0:.1:1],'Normalization','cdf');
h_dist = histcounts(hemi_sim,[0:.1:1],'Normalization','cdf');
t_dist = histcounts(temporal_sim,[0:.1:1],'Normalization','cdf');
r_dist = histcounts(temporal_unrelated,[0:.1:1],'Normalization','cdf');
plot([.05:.1:1]-.05,tc_dist,'Color',map(3,:),'LineWidth',5)
plot([.05:.1:1]-.05,h_dist,'Color',map(2,:),'LineWidth',5)
plot([.05:.1:1]-.05,t_dist,'Color',map(1,:),'LineWidth',5)
plot([.05:.1:1]-.05,r_dist,'Color','k','LineWidth',5)
xlabel('Synapse Similarity')
ylabel('Frequency')
set(gca,'FontSize',18)


%histogram(temp_cohort_sim,0:.1:1,'Normalization','probability','FaceColor',map(3,:))
%histogram(hemi_sim,0:.1:1,'Normalization','probability','FaceColor',map(1,:),'FaceAlpha',.5)
%histogram(temporal_unrelated,0:.1:1,'Normalization','probability','FaceColor','k','FaceAlpha',.1)

%%
covmap = flipud(cbrewer('seq','YlOrRd',64))
lr_cov = cov(presim_mat_l - presim_mat_r)
imagesc(lr_cov,[-.6e-3 .6e-3]); colormap(covmap)
%%
for i = 1:4
    t = tc_sim(table2array(tc_index(:,3)) == i)
    c{i} = ones(length(cat(1,t{:})),1)+i
    within_cohort_sim{i} = cat(1,t{:})
end

data = cat(1,within_cohort_sim{:})
labels = cat(1,c{:})

anova_stats = anovan(data,labels)

%% Compare NB1-2 to others

figure; hold on

barh(1,nanmean(temp_cohort_sim),'FaceColor','k','FaceAlpha',.5)
errorbar(nanmean(temp_cohort_sim),1,nanstd(temp_cohort_sim)/sqrt(numel(temp_cohort_sim)),'k','horizontal')

barh(2,nanmean(NB1_2_analysis.Synapse_Similarity.Presynapse_Similarity),'FaceColor','b','FaceAlpha',.5)
yticks([1,2])
yticklabels({'Total Cohort Similarity','NB1-2 Early Similarity'})
xlabel('Presynaptic Similarity')
set(gca,'FontSize',18)
set(gca,'YDir','reverse')
