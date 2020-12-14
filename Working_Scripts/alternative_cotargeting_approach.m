load mark2020_neurons_temporal_cohorts.mat
tc_unique = unique(an_in(:,[2,4,6]),'rows');
[~,tc_index] = ismember(an_in(:,[2,4,6]),tc_unique,'rows')

neuron_deg = get_adjacency(nl,0,2)
neuron_deg(neuron_deg>0) = 1;


for i = 1:max(tc_index)
    clear temporal_cohort_pairs temporal_cohort_pairs
    temporal_cohort_pairs = nchoosek(find(tc_index==i),2)
    if length(temporal_cohort_pairs)>1
        for j = 1:size(temporal_cohort_pairs,1)
            neurons_tc{i}{j} = nl(temporal_cohort_pairs(j,:))
            
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

shared_tc_targets = cat(2,common_targeting_tc{:})
shared_tc_targets_bin = shared_tc_targets;
shared_tc_targets_bin(shared_tc_targets_bin>0) = 1
temporal_cohort_cotarget_probability = sum(shared_tc_targets_bin)/numel(shared_tc_targets_bin)
%% 
sim_thresh = [0:.05:.2]
for s_counter = 1:length(sim_thresh)
    sim_min = sim_thresh(s_counter)
    for kk = 1:1000
        clc
        display(kk/1000)
        tc_index_shuffle = tc_index(randperm(length(tc_index)));
        for i = 1:max(tc_index_shuffle)
            clear temporal_cohort_pairs temporal_cohort_pairs 
            if length(find(tc_index_shuffle==i))>1
                temporal_cohort_pairs = nchoosek(find(tc_index_shuffle==i),2);
                for j = 1:size(temporal_cohort_pairs,1);
                    neurons{i}{j} = nl(temporal_cohort_pairs(j,:));

                    neuron_1_connectivity = neuron_deg(temporal_cohort_pairs(j,1),:);
                    neuron_2_connectivity = neuron_deg(temporal_cohort_pairs(j,2),:);

                    neuron_1_connectivity_temporal_cohorts = tc_index(neuron_1_connectivity==1);
                    neuron_2_connectivity_temporal_cohorts = tc_index(neuron_2_connectivity==1);

                    connectivity_overlap = intersect(neuron_1_connectivity_temporal_cohorts,neuron_2_connectivity_temporal_cohorts);

                    common_targeting_random{kk}{i}(j) = numel(connectivity_overlap);
                    clear neuron_1_connectivity and neuron_2_connectivity and neuron_1_connectivity_temporal_cohorts and neuron_2_connectivity_temporal_cohorts
                end
            else
            end

        end
    end
end
%%
for i = 1:1000
    shared_tc_targets_rand = cat(2,common_targeting_random{i}{:})
    shared_tc_targets_bin_rand = shared_tc_targets_rand;
    shared_tc_targets_bin_rand(shared_tc_targets_bin_rand>0) = 1
    temporal_cohort_cotarget_probability_rand(i) = sum(shared_tc_targets_bin_rand)/numel(shared_tc_targets_bin_rand)
end

