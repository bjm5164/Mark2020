% This script will generate the connectivity statistics in for figure 4.
% First, we make an adjacency matrix from the full dataset.  Then we
% generate a combined adjacency matrix for each hemilineage, sensory, or
% motor pool.  The combined adjacency matrix weights are total outputs of
% hemilineage i onto hemilineage j divided by the total inputs for
% hemilieange j.  It uses the scripts get_adjacency and concatinate
% synapses

load mark2020_neurons.mat

% Get the adjacency matrix of all neurons
adj = get_adjacency(Neuron_List,0); 
adj(adj<2) = 0 % Threshold for single synapses
lineages = unique(an.Lineage);

% Generate a hemilineage index
clear hemilineage_index
[~,hemilineage_index] = ismember([an.DV_Index],unique([an.DV_Index],'rows'),'rows');
hemi_legend = {'Ventral','Dorsal','MN','SN'}
%% Sum_Connectivity
% For each hemilineage, or sensory/motor neuron pool, sum all of the
% connections between groups, and then divide by the total number of inputs
% for that group to get a grouped adjacency matrix with the weights being
% fractions of inputs onto each group.

comb_adj = zeros(max(hemilineage_index));
for i = 1:max(hemilineage_index);
    hli_i = find(hemilineage_index == i); % find neurons in hemilineage i
    for j = 1:max(hemilineage_index);
        hli_j = find(hemilineage_index == j); % find neurons in hemilineage j
        hemi_cat_j = concatinate_synapses(Neuron_List(hli_j),hemi_legend{j}); % Concatinate all the synapses and partner SkIDs for neurons in hemilineage j.  This will give us the total number of inputs for this hemilineage.
        comb_adj(i,j) = sum(reshape(adj(hli_i,hli_j),numel(adj(hli_i,hli_j)),1))/length(hemi_cat_j.Inputs.xyz); % Sum all of outputs of hemilineage i onto hemilineage j and divide by the number of inputs for hemilineage j.
        clear hli_j and hemi_cat_j
    end
    clear hli_i
end

figure;
output_names = transpose(hemi_legend)
input_names = transpose(hemi_legend);
map = cbrewer('seq','Blues',64);
imagesc(comb_adj,[0 .1]); xticks(1:length(input_names)); 
xticklabels({input_names{:}}); xtickangle(90); 
yticks(1:length(output_names));yticklabels({output_names{:}}); 
axis xy;
colormap(map);
set(gca,'FontSize',18);
c = colorbar;

G = digraph(comb_adj)
figure;
pg = plot(G,'Layout','force','NodeLabel',input_names,'EdgeLabel',G.Edges.Weight,'EdgeAlpha',1,'EdgeColor','k')
pg.XData = [2,2,3,1]
pg.YData = [3,1,2,2]
pg.ArrowPosition = 1
pg.ArrowSize = 200*G.Edges.Weight
pg.LineWidth = 100*G.Edges.Weight
hold on
plot([.25,.75],[.25,.25],'LineWidth',100*.05)
plot([.25,.75],[.5,.5],'LineWidth',100*.1)
plot([.25,.75],[.75,.75],'LineWidth',100*.15)
