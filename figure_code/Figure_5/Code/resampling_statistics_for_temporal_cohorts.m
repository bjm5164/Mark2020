load mark2020_neurons_temporal_cohorts.mat
%% Generate similarity matrix and grouping variables

an_sided = an_in(an_in.Side_Index == 1,:);

% Grouping vector for stem-cell
unique_hemilineage = unique(an_sided(:,[2,4]),'rows')
[~,hemilineage_grouping] = ismember(an_sided(:,[2,4]),unique_hemilineage,'rows')

% Grouping vector for stem-cell AND birth-order related
unique_temporal = unique(an_sided(:,[2,4,6]),'rows')
[~,temporal_grouping] = ismember(an_sided(:,[2,4,6]),unique_temporal,'rows')

% Calculate left/left similarity and right/right similarity and average
similarity.left = synapse_similarity(neurons(an_in.Side_Index == 0),2000,3,[],1)
similarity.right = synapse_similarity(neurons(an_in.Side_Index == 1),2000,3,[],2)
similarity.combined = (similarity.left + similarity.right)*.5
%%

% For each individual group of neurons related by stem cell and birth order
for i = 1:max(temporal_grouping)
    related_index = temporal_grouping == i % Get index of related neurons
    unrelated_index = boolean((related_index*-1)+1) % Get index of unrelated neurons
    
    sim_mat_related = similarity_matrix(related_index,related_index) % Get subset of related comparison
    sim_mat_unrelated = similarity_matrix(related_index,unrelated_index) % Get subset of unrelated comparisons
    
    related_vals = sim_mat_related(boolean(triu(ones(size(sim_mat_related)),1))) % Similarity matrix is symmetric, so get each pairwise comparison only once.
    unrelated_vals = sim_mat_unrelated(:) % Unrelated matrix isn't, so just concatinate all values.
    
    related{i} = related_vals % Store in cell since each group is not the same size
    unrelated{i} = unrelated_vals
end


intra_cohort = cat(1,related{:});
inter_cohort = cat(1,unrelated{:}); clear related and unrelated

% Calculate real difference
real_mean_dif = mean(intra_cohort) - mean(inter_cohort) 

%%
for j = 1:10000
   
    % Initialize indexing
    an_rand = an_sided;
    % For each hemilineage, randomize the birth order
    unique_hemilineage = unique(an_sided(:,[2,4]),'rows');
    for i = 1:height(unique_hemilineage);
        hemi_index = find(ismember(an_sided(:,[2,4]),unique_hemilineage(i,:),'rows'));
        tc = an_sided(hemi_index,6);
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