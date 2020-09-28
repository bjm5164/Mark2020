test_skids = neurons_used_rand{4}

neuron_pairs = reshape(test_skids,1000*100,2);
unique_pairs = unique(neuron_pairs,'rows')

n1 = find([nl(:).SkIDs] == unique_pairs(end,1))
n2 = find([nl(:).SkIDs] == unique_pairs(end,2))

figure; hold on
for i = 1:length(similarity_distributions_random)
    bar(i,mean(similarity_distributions_random{i}(:)),'FaceColor',map(i,:))
end

bar(i+1,mean(similarity_distributions(:)),'FaceColor','k')


%%
sim_mat = synapse_similarity_v2(nl,2000,3,[],1)

full_adj = get_adjacency(Neuron_List,1)
full_adj(full_adj<.01) = 0;
full_adj(an.DV_Index>1,:) = []

connectivity_distance = 1-squareform(pdist(adj_real,'cosine'))
figure;
scatter(sim_mat(:),connectivity_distance(:))
%%
adj_to_use = adj_frac

for i = 1:length(adj_to_use)
    for j = 1:length(adj_to_use)
        con_sim(i,j) = adj_to_use(i,j)/sum(adj_to_use(i,:))*(adj_to_use(i,j)/sum(adj_to_use(:,j)));
    end
end

figure;
scatter(sim_mat(:),con_sim(:))

%%
figure; hold on
histogram(co_target_fraction(:),0:.1:1,'FaceColor','k')
%ctf_real = histcounts(mean(co_target_fraction),0:.1:1)
%area(0:.1:.9,ctf_real)
for i = 1:length(co_target_fraction_rand)
    histogram(co_target_fraction_rand{i}(:),0:.1:1,'FaceColor',map(i,:))
    %ctf_rand = histcounts(mean(co_target_fraction_rand{i}),0:.1:1)
    %area(0:.1:.9,ctf_rand)
end

%%
figure; hold on
plot_neurons(nl(111),'k',1,3,1,1,0,[],{'r','c'})
plot_neurons(nl(115),'k',1,3,1,1,0,[],{'m','b'})

