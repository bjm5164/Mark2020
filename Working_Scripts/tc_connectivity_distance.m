adj(boolean(eye(size(adj)))) = 0
[comb_adj,input_size] = combined_adj(adj>0,tc_index)

%%


real_intra = comb_adj(boolean(eye(size(comb_adj))))

connect_strengths_real = mean(comb_adj(comb_adj>0));
density_real = sum(comb_bin(:)) ./ numel(comb_bin);

%%

for i = 1:length(sim_thresh)
    for j = 1:1000
        adj_shuf = Wshuf_sim(:,:,j,i);
        comb_shuf = combined_adj(adj_shuf,tc_index);
        shuf_intra = comb_shuf(boolean(eye(size(comb_shuf))));
        
        intra_prob_shuf(:,j,i) = shuf_intra;
        connect_strengths_shuf(i,j) = mean(comb_shuf(comb_shuf>0));
        density_shuf(i,j) = sum(comb_shuf(:)) ./ numel(comb_shuf);
        
        if i == 1
            adj_shuf_deg = Wshuf_deg(:,:,j,i);
            comb_shuf_deg = combined_adj(adj_shuf_deg,tc_index);
            density_deg_shuf(j) = sum(comb_shuf_deg(:)) ./ numel(comb_shuf_deg);
        end
    end
end
        
%% Intra-cohort connectivity probability
cmap = cool(3);
figure; hold on
for i = 1:3
    counts = histcounts(mean(intra_prob_shuf(:,:,i)),[0:.005:.05],'Normalization','probability')
    area([.0025:.005:.05],counts,'FaceColor',cmap(:,i),'FaceAlpha',.5)
end
mean(real_intra)
std(real_intra)

%% Connection Strength (weights of non-zero cohort connections)
figure; hold on
for i = 1:3
    histogram(connect_strengths_shuf(i,:),[.2:.0025:.35],'Normalization','probability','FaceColor',cmap(:,i))
end

%% Connection Density
figure; hold on
for i = 1:3
    histogram(density_shuf(i,:),'FaceColor',cmap(:,i))
end


%%
function [comb_adj,input_size] = combined_adj(adj,index)
   
    comb_adj = zeros(max(index));
    input_size = zeros(max(index));
    for i = 1:max(index);
        group_ind_i = find(index == i);
        for j = 1:max(index);
            group_ind_j = find(index == j);
            comb_adj(i,j) = sum(reshape(adj(group_ind_i,group_ind_j),numel(adj(group_ind_i,group_ind_j)),1));
            input_size(i,j) = numel(group_ind_j);
            clear group_ind_j
        end
        clear group_ind_i
    end
end