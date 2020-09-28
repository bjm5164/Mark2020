load mark2020_neurons_temporal_cohorts.mat
adj_frac = get_adjacency(nl,1);
adj_frac(adj_frac<.01) = 0;
neuron_deg = adj_frac;
neuron_deg(neuron_deg>0) = 1;
sim_mat = synapse_similarity_v2(nl,2000,3,[],1);
sim_mat(boolean(eye(size(sim_mat)))) = 0;
tc_unique = unique(an_in(:,[2,4,5,6]),'rows');
[~,tc_index] = ismember(an_in(:,[2,4,5,6]),tc_unique,'rows')

%% First, look at the distribution of temporal cohort sizes.
[C,ia,ic] = unique(tc_index);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];
figure; histogram(a_counts,'Normalization','probability'); hold on

% Resample to generate a distribution of cohort sizes
%rdist = a_counts(randsample(length(a_counts),1000,'true'))
%histogram(rdist,'Normalization','probability')
%%


figure; hold on
histogram(tc_index,[0:1:75],'Normalization','probability');
rdist = tc_index(randsample(length(tc_index),1000,'true'))
histogram(rdist,[0:1:75],'Normalization','probability');
%%
for i = 1:max(tc_index)
    clear temporal_cohort_pairs temporal_cohort_pairs
    temporal_cohort_pairs = nchoosek(find(tc_index==i),2)
    tcp{i} = temporal_cohort_pairs
    if length(temporal_cohort_pairs)>1
        for j = 1:size(temporal_cohort_pairs,1)
            %neurons_tc{i}{j} = nl(temporal_cohort_pairs(j,:))
            tc_similarity{i}(j) = sim_mat(temporal_cohort_pairs(j,1),temporal_cohort_pairs(j,2));
            neuron_1_connectivity = neuron_deg(temporal_cohort_pairs(j,1),:)
            neuron_2_connectivity = neuron_deg(temporal_cohort_pairs(j,2),:)
            
            neuron_1_connectivity_temporal_cohorts = tc_index(neuron_1_connectivity==1);
            neuron_2_connectivity_temporal_cohorts = tc_index(neuron_2_connectivity==1);
            
            connectivity_overlap = intersect(neuron_1_connectivity_temporal_cohorts,neuron_2_connectivity_temporal_cohorts);
            
            common_targeting_tc{i}(j) = numel(connectivity_overlap)
            
            clear neuron_1_connectivity and neuron_2_connectivity and neuron_1_connectivity_temporal_cohorts and neuron_2_connectivity_temporal_cohorts
        end
    else
    end
   
end

shared_targets_tc = cat(2,common_targeting_tc{:})
shared_targets_tc_bin = shared_targets_tc;
shared_targets_tc_bin(shared_targets_tc_bin>0) = 1
cotarget_probability_temporal_cohort = sum(shared_targets_tc_bin)/numel(shared_targets_tc_bin)

%%

sim_thresh = [.01:.05:.3];
for i = 1:length(sim_thresh)
    [row,col] = ind2sub(size(sim_mat),find(sim_mat>sim_thresh(i)));
    sample_pairs = [row,col];

    for ii = 1:length(sample_pairs)
         if tc_index(sample_pairs(ii,1)) == tc_index(sample_pairs(ii,2))
            common_targeting_prox{i}(ii) = nan;
         else
            neurons_prox{i}{ii} = nl([sample_pairs(ii,1),sample_pairs(ii,2)]);

            neuron_1_connectivity = neuron_deg(sample_pairs(ii,1),:);
            neuron_2_connectivity = neuron_deg(sample_pairs(ii,2),:);

            neuron_1_connectivity_temporal_cohorts = tc_index(neuron_1_connectivity==1);
            neuron_2_connectivity_temporal_cohorts = tc_index(neuron_2_connectivity==1);

            connectivity_overlap = intersect(neuron_1_connectivity_temporal_cohorts,neuron_2_connectivity_temporal_cohorts);

            common_targeting_prox{i}(ii) = numel(connectivity_overlap);

         end
        clear neuron_1_connectivity and neuron_2_connectivity and neuron_1_connectivity_temporal_cohorts and neuron_2_connectivity_temporal_cohorts
    end
end


%%
figure; hold on
histogram(shared_targets_tc,'Normalization','cdf','FaceColor','m')
colors = cbrewer('seq','Blues',4)
for i = 1:length(common_targeting_prox)
    shared_targets_prox = common_targeting_prox{i};
    shared_targets_prox_bin = shared_targets_prox
    shared_targets_prox_bin(shared_targets_prox_bin>0) = 1;
    cotarg_sums = nansum(shared_targets_prox_bin,2)
    

    cotarget_probability_prox(i) = nansum(shared_targets_prox_bin)/sum(isnan(shared_targets_prox_bin) == 0)
end
