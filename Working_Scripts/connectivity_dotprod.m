
real_co_target_dist = cohort_target_dist(adj,tc_index)

for i = 1:3
    shuf_cohort_target_dist(:,:,i) = cohort_target_dist(Wshuf_sim(:,:,:,i),tc_index)
end

shuf_cohort_target_dist_deg = cohort_target_dist(Wshuf_deg(:,:,:,1),tc_index)


%%
figure; hold on
for i = 1:3
    counts = histcounts(mean(shuf_cohort_target_dist(:,:,i),1),[0:.02:5],'Normalization','probability')
    area([.01:.02:5],counts,'FaceColor',map(i,:),'FaceAlpha',.5); 
end

counts = histcounts(mean(shuf_cohort_target_dist_deg,1),[0:.02:5],'Normalization','probability')
area([.01:.02:5],counts,'FaceColor','k','FaceAlpha',.5);

plot([mean(real_co_target_dist), mean(real_co_target_dist)],[0, 1],'Color','r','LineWidth',3,'LineStyle','--')
xlabel('Mean Dot Product')
ylabel('Frequency')
set(gca,'FontSize',18)


%%
function tc_connectivity_dist = cohort_target_dist(adj_set,index)
iter_number = size(adj_set,3);

for iter = 1:iter_number
    adj_iter = adj_set(:,:,iter);
    adj_iter = adj_iter';
    
    for i = 1:max(index)
        clear temporal_cohort_pairs temporal_cohort_pairs
        temporal_cohort_pairs = nchoosek(find(index==i),2);
        if length(temporal_cohort_pairs)>1
            for j = 1:size(temporal_cohort_pairs,1);

                neuron_1_connectivity = adj_iter(temporal_cohort_pairs(j,1),:);
                neuron_2_connectivity = adj_iter(temporal_cohort_pairs(j,2),:);

                dist{i}(j) = dot(neuron_1_connectivity,neuron_2_connectivity);
                %dist{i}(j) = pdist2(neuron_1_connectivity,neuron_2_connectivity,'correlation');

             
            end
        else
        end

    end

    target_distances = cat(2,dist{:});
    target_distances(isnan(target_distances)) = 1;

    tc_connectivity_dist(:,iter) = target_distances; clear target_distances
end
end