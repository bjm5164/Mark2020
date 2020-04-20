load m19_neurons.mat

adj = get_adjacency(nl,1)

tc_unique = unique(an_in(:,[2,4,6]),'rows');
[~,tc_index] = ismember(an_in(:,[2,4,6]),tc_unique,'rows')

temporal_cohorts = arrayfun(@(x) nl(find(tc_index==x)),1:max(tc_index),'UniformOutput',false)

for i = 1:length(temporal_cohorts)
    tc_cat(i) = concatinate_synapses(temporal_cohorts{i},strcat(unique(an_in.Lineage(tc_index==i)),'_Cohort:',num2str(unique(an_in.Temporal_Cohort(tc_index==i)))))
end

%%

[sim_mat] = synapse_similarity_io_overlap(nl,2000)

neuron_deg = get_adjacency(nl,1);
neuron_deg(neuron_deg<.01) = 0;
%neuron_deg(neuron_deg>0) = 1;
neuron_adj_ud = neuron_deg+neuron_deg';


%%
for kk = 1:1000
    for i = 1:100
        sim = 0;
        while sim<.2
            % Choose a temporal cohort.  Make sure the temporal
            % cohort has at least two neurons in it. Then pick a pair of
            % neurons and get their connectivity.  Make sure the first neuron
            % has at least some connectivity to other neurons in the dataset.
            neuron_ip_1_connections = 0;
            while neuron_ip_1_connections<1 
                cohort_ip_index = 0;
                while length(cohort_ip_index) < 2
                    cohort_ip = randperm(max(tc_index),1);
                    cohort_ip_index = find(tc_index==cohort_ip);     
                end
                neurons_ip = cohort_ip_index(randperm(length(cohort_ip_index),2));
                % Get the connectivity of neuron_ip_1 and make sure it has at least one
                % connection in our dataset.
                neuron_ip_1_connection_index = find(neuron_deg(neurons_ip(1),:) > 0);
                neuron_ip_1_connections = length(neuron_ip_1_connection_index);
            end
            sim = sim_mat(neurons_ip(1),neurons_ip(2));
        end

        % Pick a neuron that is connected to neuron ip, we'll call it
        % neuron_jq.  Get the temporal cohort of neuron_jq.
        neuron_jq = neuron_ip_1_connection_index(randperm(length(neuron_ip_1_connection_index),1));
        cohort_jq = tc_index(neuron_jq);

        % Get connectivity of neuron_ip_2
        neuron_ip_2_connection_index = find(neuron_deg(neurons_ip(2),:) > 0);
        neuron_ip_2_connection_cohorts = tc_index(neuron_ip_2_connection_index);

        % Check if any of neuron_ip_2 connections are also of cohort_jq
        same_cohort(kk,i) = ismember(cohort_jq,neuron_ip_2_connection_cohorts);

    end
    
end  


for kk = 1:1000
    display(kk)
    for i = 1:100
        sim = 0;
        while sim<.2
            tc_index_randomized = tc_index(randperm(length(tc_index)));
            % Choose a temporal cohort.  Make sure the temporal
            % cohort has at least two neurons in it. Then pick a pair of
            % neurons and get their connectivity.  Make sure the first neuron
            % has at least some connectivity to other neurons in the dataset.
            neuron_ip_1_connections = 0;
            while neuron_ip_1_connections<1
                cohort_ip_index = 0;
                while length(cohort_ip_index) < 2
                    cohort_ip = randperm(max(tc_index_randomized),1);
                    cohort_ip_index = find(tc_index_randomized==cohort_ip);
                    if length(unique(tc_index(cohort_ip_index)))<2
                        cohort_ip_index = 0;
                    else
                    end
                end
                neurons_ip = cohort_ip_index(randperm(length(cohort_ip_index),2));
                % Get the connectivity of neuron_ip_1 and make sure it has at least one
                % connection in our dataset.
                neuron_ip_1_connection_index = find(neuron_deg(neurons_ip(1),:) > 0);
                neuron_ip_1_connections = length(neuron_ip_1_connection_index);
            end
            sim = sim_mat(neurons_ip(1),neurons_ip(2));
        end
        % Pick a neuron that is connected to neuron ip, we'll call it
        % neuron_jq.  Get the temporal cohort of neuron_jq.
        neuron_jq = neuron_ip_1_connection_index(randperm(length(neuron_ip_1_connection_index),1));
        cohort_jq = tc_index(neuron_jq);

        % Get connectivity of neuron_ip_2
        neuron_ip_2_connection_index = find(neuron_deg(neurons_ip(2),:) > 0);
        neuron_ip_2_connection_cohorts = tc_index(neuron_ip_2_connection_index);

        % Check if any of neuron_ip_2 connections are also of cohort_jq
        same_cohort_random(kk,i) = ismember(cohort_jq,neuron_ip_2_connection_cohorts);

    end
    
end  
sum(mean(same_cohort,2))/10000
sum(mean(same_cohort_random,2))/10000

figure; hold on
histogram(mean(same_cohort,2),0:.02:1)
histogram(mean(same_cohort_random,2),0:.02:1)
%% Sum_Connectivity
comb_adj = zeros(length(tc_index));
for i = 1:length(tc_index)
    tc_i = tc_index{i};
    for j = 1:length(tc_index)
        tc_j = tc_index{j};
        comb_adj(i,j) = sum(reshape(neuron_deg(tc_i,tc_j),numel(neuron_deg(tc_i,tc_j)),1));
        clear tc_j 
    end
    clear tc_i
end

%%
map = hsv(max(tc_unique.Lineage_Index))
for i = 30:33
figure
gscatter(sim_mat(i,:),comb_adj(i,:),tc_unique.Lineage_Index(:))
title(tc_cat(i).Name)
xlabel('similarity')
ylabel('connectivity')
end


%%
map = hsv(6)
 
Synapse_Distributions(nl(1),'r',1,.25,0)
for i = 1:6 
    Synapse_Distributions(nl(i),map(i,:),2,.25,0)
end
colorbar; colormap(map)
    

[ii,jj] = ind2sub(size(sim_mat),find(sim_mat>.2 & comb_adj==0))

for i = 1:10
    figure; hold on
    Synapse_Distributions(tc_cat(ii(i)),'r',1,.25,0)
    Synapse_Distributions(tc_cat(jj(i)),'c',2,.25,0)
end

%%
for i = 1:length(sim_mat)
    post_cohort_targets(i) = sum(sim_mat(:,i)>.3);
end
figure
histogram(post_cohort_targets)
xlabel('Cohorts Targeting a Single Cohort')

load m19_neurons.mat
[sim_mat_all] = synapse_similarity_io_overlap(nl,1000)
adj = get_adjacency(nl,1)


%%
tc_unique = unique(an_in(:,[5,2,4,6]),'rows')
tc_index = arrayfun(@(x) ismember(an_in(:,[5,2,4,6]),tc_unique(x,:),'rows'),1:size(tc_unique,1),'UniformOutput',false)

for side = 0:1
    for dv = 0:1
        for lin = 1:max(an_in.Lineage_Index(an_in.DV_Index == dv & an_in.Side_Index == side))
            for tc = 1:max(an_in.Temporal_Cohort(an_in.Side_Index == side & an_in.DV_Index == dv & an_in.Lineage_Index == lin))
                tc_adj = adj(an_in.Side_Index == side & an_in.DV_Index == dv & an_in.Lineage_Index == lin & an_in.Temporal_Cohort == tc,:)
                for kk = 1:size(tc_adj,1)
                    lineage_relationship{side+1}{dv+1}{lin}{tc}{kk} = [an_in.Side_Index(tc_adj(kk,:) > .01), an_in.DV_Index(tc_adj(kk,:) > .01),an_in.Lineage_Index(tc_adj(kk,:)>.01),an_in.Temporal_Cohort(tc_adj(kk,:)>.01)]
                end
            end
        end
    end
end
