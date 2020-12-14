load mark2020_neurons_temporal_cohorts.mat
%%
an_half = an_in(1:2:end,:)

lineage_grouping = an_in.Lineage_Index(an_in.Side_Index == 0)

unique_hemilineage = unique(an_half(:,[2,4]),'rows')
[~,hemilineage_grouping] = ismember(an_half(:,[2,4]),unique_hemilineage,'rows')

unique_temporal = unique(an_half(:,[2,4,6]),'rows')
[~,temporal_grouping] = ismember(an_half(:,[2,4,6]),unique_temporal,'rows')


similarity.left = synapse_similarity_v2(nl(an_in.Side_Index == 0),2000,3,[],1)
similarity.right = synapse_similarity_v2(nl(an_in.Side_Index == 1),2000,3,[],2)
similarity.combined = (similarity.left + similarity.right)*.5
%%
[~,temporal_grouping] = ismember(an_half(:,[2,4,6]),unique_temporal,'rows')

for i = 1:max(temporal_grouping)
    related_index = temporal_grouping == i
    unrelated_index = boolean((related_index*-1)+1)
    
    sim_mat_related = similarity.combined(related_index,related_index)
    sim_mat_unrelated = similarity.combined(related_index,unrelated_index)
    
    related_vals = sim_mat_related(boolean(triu(ones(size(sim_mat_related)),1)))
    unrelated_vals = sim_mat_unrelated(:)
    
    related{i} = related_vals
    unrelated{i} = unrelated_vals
end


intra_lineage = cat(1,related{:});
inter_lineage = cat(1,unrelated{:}); clear related and unrelated

real_mean_dif = mean(intra_lineage) - mean(inter_lineage)

%%
for j = 1:10000
   
    % Initialize indexing
    an_rand = an_half;
    % For each hemilineage, randomize the birth order
    unique_hemilineage = unique(an_half(:,[2,4]),'rows');
    for i = 1:height(unique_hemilineage);
        hemi_index = find(ismember(an_half(:,[2,4]),unique_hemilineage(i,:),'rows'));
        tc = an_half(hemi_index,6);
        an_rand(hemi_index,6) = tc(randperm(height(tc)),:); clear tc and hemi_index
    end
    
    % Get index for unique birth cohorts
    unique_temporal = unique(an_rand(:,[2,4,6]),'rows');
    [~,temporal_grouping] = ismember(an_rand(:,[2,4,6]),unique_temporal,'rows');
    
    for i = 1:max(temporal_grouping)
        related_index = temporal_grouping == i;
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
    
    random_mean_dif(j) = mean(intra_random) - mean(inter_random);
    
end
    
figure; hold on
histogram(random_mean_dif)

%%
stats.lineage = f_anosim(1-similarity.combined,lineage_grouping,[],10000,0,1)


stats.temporal = f_anosim(1-similarity.combined,temporal_grouping,[],10000,0,1)
%%

stats.hemilineage = f_anosim(1-similarity.combined,hemilineage_grouping,[],10000,0,1)

%% dv grouping
[~,dv_grouping] = ismember(an_half.DV_Index,[0,1])
stats.dv = f_anosim(1-similarity.combined,dv_grouping,[],10000,0,1)

%%
clear temporal_cohort_index and hl_mat

for i = 1:max(hemilineage_grouping)
    hemilineage = i

    temporal_cohort_index = temporal_grouping(hemilineage_grouping == hemilineage)
    hl_mat = similarity.combined(hemilineage_grouping == hemilineage,hemilineage_grouping == hemilineage)

    result.p = f_anosim(1-hl_mat,temporal_cohort_index,0,10000,0,1)
    p(i) = result.p
end