load mark2020_neurons_temporal_cohorts.mat
hemilineage_index = unique([an_in.Lineage_Index,an_in.DV_Index],'rows')

hemilineages = arrayfun(@(x) nl(an_in.Lineage_Index == hemilineage_index(x,1) & an_in.DV_Index == hemilineage_index(x,2)), 1:length(hemilineage_index),'UniformOutput',false)

for i = 1:length(hemilineage_index)
    if hemilineage_index(i,2) == 0
      hemi_legend{i} = strcat(unique(an_in.Lineage(an_in.Lineage_Index == hemilineage_index(i,1))),'Ventral')
    elseif hemilineage_index(i,2) == 1
      hemi_legend{i} = strcat(unique(an_in.Lineage(an_in.Lineage_Index == hemilineage_index(i,1))),'Dorsal')
    else
      hemi_legend{i} = unique(an_in.Lineage(an_in.Lineage_Index == hemilineage_index(i,1)))
    end
end 

dorsal_hemilineages = hemilineages(hemilineage_index(:,2) == 1);
ventral_hemilineages = hemilineages(hemilineage_index(:,2) == 0);
clearvars -except dorsal_hemilineages and ventral_hemilineages and hemilineages and Neuron_List and savefigs and directory an an and hemi_legend and nl and an_in

%% Quantify similarity

presim_mat_l = synapse_similarity_v2(nl(an_in.Side_Index == 0),2000,3,[],1);
presim_mat_r = synapse_similarity_v2(nl(an_in.Side_Index == 1),2000,3,[],1);
presim = .5*(presim_mat_l+presim_mat_r)

postsim_mat_l = synapse_similarity_v2(nl(an_in.Side_Index == 0),2000,3,[],2);
postsim_mat_r = synapse_similarity_v2(nl(an_in.Side_Index == 1),2000,3,[],2);
postsim = .5*(postsim_mat_l+postsim_mat_r)

sim_mat = (presim+postsim)*.5
d_edge_total = [nl(an_in.Side_Index == 0).Temporal_Cohort];
an_half = an_in(1:2:end,:);
nl_half_skids = [nl(1:2:end).SkIDs];

%% Get the similarity scores for all neurons related by birthdate, hemilineage, or temporal cohort

% Define which similarity matrix to use
similarity_matrix = presim

% Temporal
temporal_index = unique(an_half(:,6)); % Possible temporal cohorts

for i = 1:length(temporal_index.Temporal_Cohort)
    related_index = ismember(an_half(:,6),temporal_index(i,:)) % Find neurons in temporal cohort i
    sr = similarity_matrix(related_index,related_index) %get the pairwise similarities between them
    sur = similarity_matrix(related_index,boolean(related_index*-1+1)) % Get pairwise comparisons of neurons in temporal cohort i and all neurons not in temporal cohort i
    sim_related = sr(:)
    sim_unrelated = sur(:)
    sim_related(sim_related == 1) = [] % remove self-self comparisons
    t_sim{i} = sim_related; 
    nt_sim{i} = sim_unrelated; clear sim_related and sim_unrelated and sr and sur
end

% Hemilineage 
hl_index = unique(an_half(:,[2,4]),'rows') % Possible hemilineages

for i = 1:length(hl_index.DV_Index)
    related_index_h = ismember(an_half(:,[2,4]),hl_index(i,:),'rows'); % Find neurons in hemilineage i
    sr_h = similarity_matrix(related_index_h,related_index_h); % get pairwise similarities between them
    sur_h = similarity_matrix(related_index_h,boolean((related_index_h*-1)+1)); % Get pairwise comparisons of neurons in hemilineage i and all neurons not in hemilineage i
    sim_related = sr_h(:);
    sim_unrelated = sur_h(:);
    sim_related(sim_related == 1) = []; % remove self-self comparisons
    h_sim{i} = sim_related; 
    nh_sim{i} = sim_unrelated; clear sim_related and sim_unrelated and sr and sur
end

% Temporal Cohort
tc_index = unique(an_half(:,[2,4,6]),'rows') % Possible temporal cohorts

for i = 1:length(tc_index.DV_Index)
    related_index = ismember(an_half(:,[2,4,6]),tc_index(i,:),'rows') % Find neurons in temporal cohort i
    sr = similarity_matrix(related_index,related_index)
    sur = similarity_matrix(related_index,boolean((related_index*-1)+1));
    sim_unrelated = sur(:);
    sim_related = sr(:);
    sim_related(sim_related == 1) = []
    tc_sim{i} = sim_related;
    ntc_sim{i} = sim_unrelated; clear sim_related and sim_unrelated and sr and sur
end


for i = 1:length(tc_sim)
    if isempty(tc_sim{i})
        idx(i) = 1
    else
        idx(i) = 0
    end
end

tc_sim(boolean(idx)) = []
%%
temporal_sim = cat(1,t_sim{:});
temporal_unrelated = cat(1,nt_sim{:});

hemi_sim = cat(1,h_sim{:});
hemi_unrelated = cat(1,nh_sim{:});

temp_cohort_sim = cat(1,tc_sim{:});
temp_cohort_unrelated = cat(1,ntc_sim{:});

temporal_sim(isnan(temporal_sim)) = []
hemi_sim(isnan(hemi_sim)) = []
temp_cohort_sim(isnan(temp_cohort_sim)) = []

[p_temporal_vs_nontemporal,h,stat] = ranksum(temporal_sim,temporal_unrelated)
[p_hemilineage_vs_temporal,h,stat] = ranksum(temporal_sim,hemi_sim)
[p_hemilineage_vs_temporal_cohort,h,stat] = ranksum(hemi_sim,temp_cohort_sim)
nanmean(temp_cohort_sim)
nanmean(hemi_sim)
nanmean(temporal_sim)

figure; hold on
%barh(1,nanmean(temporal_unrelated),'FaceColor',[0 0 0])
%errorbar(nanmean(temporal_unrelated),1,nanstd(temporal_unrelated)/sqrt(numel(temporal_unrelated)),'k','horizontal')
Violin(temporal_unrelated,1,'ViolinColor',[0,0,0],'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])

%barh(2,nanmean(temporal_sim),'FaceColor',[.25 .25 .25])
%errorbar(nanmean(temporal_sim),2,nanstd(temporal_sim)/sqrt(numel(temporal_sim)),'k','horizontal')
Violin(temporal_sim,2,'ViolinColor',map(1,:),'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])

%barh(3,nanmean(hemi_sim),'FaceColor',[.5 .5 .5])
%errorbar(nanmean(hemi_sim),3,nanstd(hemi_sim)/sqrt(numel(hemi_sim)),'k','horizontal')
Violin(hemi_sim,3,'ViolinColor',map(2,:),'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])

%barh(4,nanmean(temp_cohort_sim),'FaceColor',[.75 .75 .75])
%errorbar(nanmean(temp_cohort_sim),4,nanstd(temp_cohort_sim)/sqrt(numel(temp_cohort_sim)),'k','horizontal')
Violin(temp_cohort_sim,4,'ViolinColor',map(3,:),'Bandwidth',.05,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])


set(gca,'YDir','reverse')
yticks([1 2 3 4])
yticklabels({'Random','Temporal','Hemilineage','HL-Temporal Cohort'})
xlabel('Synapse Similarity')
set(gca,'FontSize',18)

figure; hold on
histogram(temp_cohort_sim,0:.05:1,'Normalization','probability')
histogram(hemi_sim,0:.05:1,'Normalization','probability')
histogram(temporal_sim,0:.05:1,'Normalization','probability')



xticks([1 2 3 4])
xlim([0 5])
ylim([0 1])
xticklabels({'Random','Temporal','Hemilineage','HL-Temporal Cohort'})
ylabel('Input Connectivity Similarity')
set(gca,'FontSize',18)
