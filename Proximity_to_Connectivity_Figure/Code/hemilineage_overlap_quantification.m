load mark2020_neurons_temporal_cohorts.mat

%%
hl_unique = unique(an_in(:,[2,3,4,5,6]),'rows');
[~,hl_index] = ismember(an_in(:,[2,4,5,6]),hl_unique(:,[1,3,4,5]),'rows');

ventral_index = an_in.DV_Index == 1;
dorsal_index = an_in.DV_Index == 0;
ventral_hl = hl_index(ventral_index);
dorsal_hl = hl_index(dorsal_index);

hemilineages = arrayfun(@(x) concatinate_synapses(nl(hl_index == x),hl_unique.Lineage(x)),1:height(hl_unique))
for i = 1:length(hemilineages)
    if contains(hemilineages(i).Included_Neurons{1},'Ventral,l')
        hemilineages(i).Name = strcat(hemilineages(i).Name,'-Ventral-L')
    elseif contains(hemilineages(i).Included_Neurons{1},'Ventral,r')
        hemilineages(i).Name = strcat(hemilineages(i).Name,'-Ventral-R')
    elseif contains(hemilineages(i).Included_Neurons{1},'Dorsal,l')
        hemilineages(i).Name = strcat(hemilineages(i).Name,'-Dorsal-L')
    elseif contains(hemilineages(i).Included_Neurons{1},'Dorsal,r')
        hemilineages(i).Name = strcat(hemilineages(i).Name,'-Dorsal-R')
    else error('wrong name')
    end
    
end

sim_mat = synapse_similarity_io_overlap(hemilineages,2000,2000)

%%

sim_mat(boolean(eye(size(sim_mat)))) = 1
Synapse_Distance_Clustering_v2(sim_mat,[hemilineages(:).Name],0,0)
%%
figure; hold on
histogram(sim_mat,[0:.1:1],'FaceColor','r')
%%
for i = 1:100
   %Initialize randomized index
    hl_index_randomized = zeros(length(hl_index),1)
    %Generate randomization for dorsal hemilineages
    random_ind = randperm(length(ventral_hl))
    hl_index_randomized(ventral_index) = ventral_hl(random_ind); clear random_ind
    %Generate randomization for ventral hemilineages
    random_ind = randperm(length(dorsal_hl))
    hl_index_randomized(dorsal_index) = dorsal_hl(random_ind); clear random_ind
      
    fake_hemilineages{i} = arrayfun(@(x) concatinate_synapses(nl(hl_index_randomized == x),hl_unique.Lineage(x)),1:height(hl_unique))
    sim_mat_random(:,:,i) = synapse_similarity_io_overlap(fake_hemilineages{i},2000,2000)
end
%%
for i = 1:size(sim_mat_random,3)
    iter = sim_mat_random(:,:,i);
    non_self = iter(~boolean(eye(size(iter))));
    all_counts_random(:,i) = histcounts(non_self,[0:.1:1],'Normalization','cdf'); clear iter and non_self
end

all_counts_real = histcounts(sim_mat(~boolean(eye(size(sim_mat)))),[0:.1:1],'Normalization','cdf')


%%

for i = 1:size(sim_mat_random,3)
    iter = sim_mat_random(:,:,i);
    iter(boolean(eye(size(iter)))) = -1;
    max_val_random(:,i) = max(iter,[],2); clear iter
end

mean_max_vals = mean(max_val_random)

for i = 1:size(max_val_random,2)
    max_counts_rand(:,i) = histcounts(max_val_random(:,i),[0:.1:1],'Normalization','cdf')
end

sim_mat_non_self = sim_mat
sim_mat_non_self(boolean(eye(size(sim_mat)))) = 0
max_counts_real = histcounts(max(sim_mat_non_self,[],2),[0:.1:1],'Normalization','cdf')

%%
figure; hold on
subplot(1,2,1)
plot([.05:.1:1],all_counts_real,'Color','r')
shadedErrorBar([.05:.1:1],mean(all_counts_random,2),std(all_counts_random,[],2))
xlabel('Pre-Post Similarity')
ylabel('Frequency')
title('Mean Hemilineage Overlap')
set(gca,'FontSize',18)
subplot(1,2,2)
plot([.05:.1:1],max_counts_real,'Color','r')
shadedErrorBar([.05:.1:1],mean(max_counts_rand,2),std(max_counts_rand,[],2))
xlabel('Pre-Post Similarity')
ylabel('Frequency')
title('Max Hemilineage Overlap')
set(gca,'FontSize',18)


%%
for i = 1:length(hemilineages)-1
    figure; hold on
    Synapse_Distributions(hemilineages(23),'r',1,.25,0)
    Synapse_Distributions(hemilineages(23),'c',2,.25,0)
end