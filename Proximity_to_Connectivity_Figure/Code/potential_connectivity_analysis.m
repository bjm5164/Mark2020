tc_unique = unique(an_in(:,[2,4,6]),'rows');
[~,tc_index] = ismember(an_in(:,[2,4,6]),tc_unique,'rows')

adj = get_adjacency(nl,0)
neuron_deg = adj;
neuron_deg(neuron_deg<2) = 0;
neuron_deg(neuron_deg>0) = 1;

potential_deg = potential_adjacency;
potential_deg(potential_deg<2) = 0;
potential_deg(potential_deg>0) = 1;

%% Real Adjacency
for kk = 1:1000
    display(kk)
    for i = 1:100
        deg = -1;
        while deg<0
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
                neuron_ip_1_connection_cohorts = tc_index(neuron_ip_1_connection_index);
            end
            deg = neuron_deg(neurons_ip(1),neurons_ip(2));
        end

        % Get connectivity of neuron_ip_2
        neuron_ip_2_connection_index = find(neuron_deg(neurons_ip(2),:) > 0);
        neuron_ip_2_connection_cohorts = tc_index(neuron_ip_2_connection_index);

        % Check if any of neuron_ip_2 connections are also of cohort_jq
        cohort_intersection = intersect(neuron_ip_1_connection_cohorts,neuron_ip_2_connection_cohorts);
        if length(cohort_intersection) > 0
            same_cohort(kk,i) = 1;
        else
            same_cohort(kk,i) = 0;
        end
        

    end
    
end  
clearvars -except same_cohort and tc_index and nl and neuron_deg and potential_deg and an_in and potential_adjacency

%% Potential adjacency
for kk = 1:1000
    display(kk)
    for i = 1:100
        deg = -1;
        while deg<0
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
                neuron_ip_1_connection_index = find(potential_deg(neurons_ip(1),:) > 0);
                neuron_ip_1_connections = length(neuron_ip_1_connection_index);
                neuron_ip_1_connection_cohorts = tc_index(neuron_ip_1_connection_index);
            end
            deg = potential_deg(neurons_ip(1),neurons_ip(2));
        end

        % Get connectivity of neuron_ip_2
        neuron_ip_2_connection_index = find(potential_deg(neurons_ip(2),:) > 0);
        neuron_ip_2_connection_cohorts = tc_index(neuron_ip_2_connection_index);

        % Check if any of neuron_ip_2 connections are also of cohort_jq
        cohort_intersection = intersect(neuron_ip_1_connection_cohorts,neuron_ip_2_connection_cohorts);
        if length(cohort_intersection) > 0
            potential_cohort(kk,i) = 1;
        else
            potential_cohort(kk,i) = 0;
        end
        

    end
    
end  
clearvars -except same_cohort and tc_index and nl and neuron_deg and potential_deg and an_in and potential_adjacency amd potential_cohort

%% Potential adjacency, random cohorts
for kk = 1:1000
        display(kk)
        for i = 1:100
            deg = -1;
            while deg<0
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
                    neuron_ip_1_connection_index = find(potential_deg(neurons_ip(1),:) > 0);
                    neuron_ip_1_connections = length(neuron_ip_1_connection_index);
                    neuron_ip_1_connection_cohorts = tc_index(neuron_ip_1_connection_index);
                end
                deg = potential_deg(neurons_ip(1),neurons_ip(2));
            end

            % Get connectivity of neuron_ip_2
            neuron_ip_2_connection_index = find(potential_deg(neurons_ip(2),:) > 0);
            neuron_ip_2_connection_cohorts = tc_index(neuron_ip_2_connection_index);

            % Check if any of neuron_ip_2 connections are also of cohort_jq

            cohort_intersection = intersect(neuron_ip_1_connection_cohorts,neuron_ip_2_connection_cohorts);
            if length(cohort_intersection) > 0
                rand_cohort(kk,i) = 1;
            else
                rand_cohort(kk,i) = 0;
            end
        end

    end  