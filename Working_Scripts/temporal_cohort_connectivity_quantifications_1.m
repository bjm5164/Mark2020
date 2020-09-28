load mark2020_neurons_temporal_cohorts.mat
sim_mat = synapse_similarity_v2(nl,2000,3,[],1)
neuron_deg = get_adjacency(nl,1)
neuron_deg(neuron_deg<.01) = 0
neuron_deg(neuron_deg>0) = 1
tc_unique = unique(an_in(:,[2,4,5,6]),'rows');
[~,tc_index] = ismember(an_in(:,[2,4,5,6]),tc_unique,'rows')

%%

for kk = 1:1000
    display(kk)
    for i = 1:100
        sim = -1;
        while sim<0
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
            sim = sim_mat(neurons_ip(1),neurons_ip(2));
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


sim_thresh = [0:.06:.3]
for jj = 1:length([0:.06:.3])

    for kk = 1:1000
        display(kk)
        for i = 1:100
            sim = -1;
            while sim<sim_thresh(jj)
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
                    neuron_ip_1_connection_cohorts = tc_index(neuron_ip_1_connection_index);
                end
                sim = sim_mat(neurons_ip(1),neurons_ip(2));
            end

            % Get connectivity of neuron_ip_2
            neuron_ip_2_connection_index = find(neuron_deg(neurons_ip(2),:) > 0);
            neuron_ip_2_connection_cohorts = tc_index(neuron_ip_2_connection_index);

            % Check if any of neuron_ip_2 connections are also of cohort_jq

            cohort_intersection = intersect(neuron_ip_1_connection_cohorts,neuron_ip_2_connection_cohorts);
            if length(cohort_intersection) > 0
                same_cohort_random{jj}(kk,i) = 1;
            else
                same_cohort_random{jj}(kk,i) = 0;
            end
        end

    end  
    

end
%%

map = cbrewer('seq','Blues',6)

figure; hold on
for i = 1:length(same_cohort_random)
    random_counts = histcounts(sum(same_cohort_random{i},2)/100,0:.04:1,'Normalization','probability')
    area([0:.04:.98],random_counts,'FaceColor',map(i,:),'FaceAlpha',.75,'LineWidth',1.5,'EdgeColor',map(i,:)*.5)
    %histogram(sum(same_cohort_random{i},2)/100,0:.02:1,'FaceColor',map(i,:),'Normalization','probability')
end

real_counts = histcounts(sum(same_cohort,2)/100,0:.04:1,'Normalization','probability')
plot([0:.04:.98],real_counts,'Color','k')
area([0:.04:.98],real_counts,'FaceColor','k','FaceAlpha',.75,'LineWidth',1.5,'EdgeColor','k')
%histogram(sum(same_cohort,2)/100,0:.02:1,'FaceColor','k','Normalization','probability')
c = colorbar; colormap(map)

c.Ticks = [0:.2:1]
c.TickLabels = [0:.08:.4]
c.Label.String = 'Similarity Threshold'
set(gca,'FontSize',18)
set(gcf,'Color','w')
xlabel('Mean Probability')
ylabel('Frequency')