function tc_target_prob = cohort_target_prob(adj_set,index)
iter_number = size(adj_set,3);

for iter = 1:iter_number
    adj_iter = adj_set(:,:,iter);
    
    for i = 1:max(index)
        clear temporal_cohort_pairs temporal_cohort_pairs
        temporal_cohort_pairs = nchoosek(find(index==i),2);
        if length(temporal_cohort_pairs)>1
            for j = 1:size(temporal_cohort_pairs,1);

                neuron_1_connectivity = adj_iter(temporal_cohort_pairs(j,1),:);
                neuron_2_connectivity = adj_iter(temporal_cohort_pairs(j,2),:);

                neuron_1_connectivity_temporal_cohorts = index(neuron_1_connectivity==1);
                neuron_2_connectivity_temporal_cohorts = index(neuron_2_connectivity==1);

                connectivity_overlap = intersect(neuron_1_connectivity_temporal_cohorts,neuron_2_connectivity_temporal_cohorts);

                common_targeting_tc{i}(j) = numel(connectivity_overlap);
                clear neuron_1_connectivity and neuron_2_connectivity and neuron_1_connectivity_temporal_cohorts and neuron_2_connectivity_temporal_cohorts
            end
        else
        end

    end

    shared_tc_targets = cat(2,common_targeting_tc{:});
    shared_tc_targets_bin = shared_tc_targets;
    shared_tc_targets_bin(shared_tc_targets_bin>0) = 1;
    tc_target_prob(iter) = sum(shared_tc_targets_bin)/numel(shared_tc_targets_bin);
end