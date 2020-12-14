load mark2020_neurons_temporal_cohorts.mat
an.Temporal_Cohort = zeros(length(an.DV_Index),1)
an.Temporal_Cohort(an.Lineage_Index > 1 & an.Lineage_Index < 9) = an_in.Temporal_Cohort


neuron_deg = get_adjacency(Neuron_List,1);
neuron_deg(neuron_deg<.01) = 0;
neuron_deg(neuron_deg>0) = 1;
neuron_adj_ud = neuron_deg+neuron_deg';
neuron_deg(find(sum(neuron_adj_ud) == 0),:) = [];
neuron_deg(:,find(sum(neuron_adj_ud) == 0)) = [];
Neuron_List(find(sum(neuron_adj_ud) == 0)) = [];
an(find(sum(neuron_adj_ud) == 0),:) = []
side_lookup = containers.Map([0,1],{'Left','Right'});
hemi_lookup = containers.Map([0:3],{'Ventral','Dorsal','Motor','Sensory'});
lin_lookup = containers.Map(unique(an.Lineage_Index),unique(an.Lineage));
hemi_legend = unique(an.Lineage);
lineages = unique(an.Lineage);


[~,hemilineage_index] = ismember([an.Side_Index,an.DV_Index,an.Lineage_Index],unique([an.Side_Index,an.DV_Index,an.Lineage_Index],'rows'),'rows');

savefigs = input('Save figs? 1:yes 0:no');
if savefigs == 1;
    directory = uigetdir;
else
    directory = 0;
end

%%
Gud = graph(neuron_deg+neuron_deg')
Gud.Edges.Weight = ones(length(Gud.Edges.Weight),1)
distance_matrix = distances(Gud)

%%
annotation_set = an
%annotation_set.DV_Index = an.DV_Index(randperm(height(an)))


% Get network distances between temporally related neurons
temporal_options = unique(annotation_set(:,[2,6,5]),'rows'); % Possible temporal cohorts
for i = 1:height(temporal_options)
    related_index = ismember(annotation_set(:,[2,6,5]),temporal_options(i,:));
    temp_distances = distance_matrix(related_index,related_index);
    temp_distances_ur = distance_matrix(related_index,boolean(related_index*-1+1))
    
    % Take only the single comparisons
    temp_distances = temp_distances(boolean(triu(ones(size(temp_distances)),1)))
    
    temp_distances_ur = temp_distances_ur(:); 
    temporal_distances{i} = temp_distances; clear temp_distances
    temporal_distances_ur{i} = temp_distances_ur; clear temp_distances and related_distances and temp_distances_ur
end



% Get network distances between hemilineage related neurons
hemilineage_options = unique(annotation_set(:,[2,4,5]),'rows')
% Remove sensory/motor options
hemilineage_options(hemilineage_options.DV_Index>1,:) = []

for i = 1:height(hemilineage_options)
    related_index = ismember(annotation_set(:,[2,4,5]),hemilineage_options(i,:));
    temp_distances = distance_matrix(related_index,related_index);
    temp_distances_ur = distance_matrix(related_index,boolean(related_index*-1+1));
    
    %Take only the sigle comparisons
    temp_distances = temp_distances(boolean(triu(ones(size(temp_distances)),1)))
    
    temp_distances_ur = temp_distances_ur(:); 
    hemilineage_distances{i} = temp_distances; clear temp_distances
    hemilineage_distances_ur{i} = temp_distances_ur; clear temp_distances and related_distances and temp_distances_ur
end



% Get network distances for temporal cohort related neurons
temporal_cohort_options = unique(annotation_set(:,[2,4,6,5]),'rows')
temporal_cohort_options(temporal_cohort_options.DV_Index>1,:) = []
% Remove NB1-2 from temporal cohort analysis
temporal_cohort_options(ismember(table2array(temporal_cohort_options) , [1 2 1 1], 'rows')...
                        | ismember(table2array(temporal_cohort_options) , [1 2 1 0], 'rows'),:) = []

for i = 1:height(temporal_cohort_options)
    related_index = ismember(annotation_set(:,[2,4,6,5]),temporal_cohort_options(i,:));
    temp_distances = distance_matrix(related_index,related_index);
    temp_distances_ur = distance_matrix(related_index,boolean(related_index*-1+1));
    
    %Take only the single comparisons
    temp_distances = temp_distances(boolean(triu(ones(size(temp_distances)),1)))
    
    temp_distances_ur = temp_distances_ur(:); 
    temporal_cohort_distances{i} = temp_distances; clear temp_distances
    temporal_cohort_distances_ur{i} = temp_distances_ur; clear temp_distances and related_distances and temp_distances_ur
end

% for i = 1:length(temporal_cohort_distances)
%     if isempty(temporal_cohort_distances{i})
%         idx(i) = 1
%     else
%         idx(i) = 0
%     end
% end
% 
% temporal_cohort_distances(boolean(idx)) = []; clear idx


% Get network distances for unrelated neurons. To do this, select a
% hemilineage, and take all neurons on the same side and measure the
% distance from that hemilineage to all other neurons. 

for i = 1:height(hemilineage_options)
    % Get hemilineage related index
    ingroup_index = ismember(annotation_set(:,[2,4,5]),hemilineage_options(i,:));
    % Get index of all neurons not in selected hemilineage
    outgroup_index = boolean(ingroup_index*-1+1)
    % Get index of all neurons on the same side (left/right and
    % dorsal/ventral) of selected hemilineage
    side_index = annotation_set.DV_Index == hemilineage_options.DV_Index(i) & annotation_set.Side_Index == hemilineage_options.Side_Index(i)
    
    % Index of neurons that are not in the selected hemilineage, but are on
    % the same side
    unrelated_index = find(outgroup_index + side_index == 2)  
    
    temp_distances_ur = distance_matrix(ingroup_index,unrelated_index);
    temp_distances_ur = temp_distances_ur(:); 
    unrelated_distances_ur{i} = temp_distances_ur; clear temp_distances and related_distances and temp_distances_ur
end



ur_d = cat(1,unrelated_distances_ur{:}); clear unrelated_distances_ur
temporal_d = cat(1,temporal_distances{:}); clear temporal_distances
temporal_d_ur = cat(1,temporal_distances_ur{:}); clear temporal_distances_ur
hemi_d = cat(1,hemilineage_distances{:}); clear hemilineage_distances
hemi_d_ur = cat(1,hemilineage_distances_ur{:}); clear hemilineage_distances_ur
tc_d = cat(1,temporal_cohort_distances{:}); clear temporal_cohort_distances
tc_d_ur = cat(1,temporal_cohort_distances_ur{:}); clear temporal_cohort_distances_ur
   
%%
% Plot synapses distances as a cumulative distribution
his_r = histcounts(ur_d,[.5:1:6.5],'Normalization','probability')
his_h = histcounts(hemi_d,[.5:1:6.5],'Normalization','probability')
his_tc = histcounts(tc_d,[.5:1:6.5],'Normalization','probability')
his_t = histcounts(temporal_d,[.5:1:6.5],'Normalization','probability')
his_all = histcounts(distance_matrix(:),[.5:1:6.5],'Normalization','probability')

figure; hold on
p1 = plot([1:6],cumsum(his_tc),'c','LineWidth',5)
p2 = plot([1:6],cumsum(his_h),'r','LineWidth',5)
p3 = plot([1:6],cumsum(his_t),'m','LineWidth',5)
p4 = plot([1:6],cumsum(his_r),'Color','k','LineWidth',5)
%p5 = plot([1:6],cumsum(his_all),'k','LineWidth',5)


scatter([1:6],cumsum(his_r),1000,'.','k')
scatter([1:6],cumsum(his_h),1000,'.','r')
scatter([1:6],cumsum(his_tc),1000,'.','c')
scatter([1:6],cumsum(his_t),1000,'.','m')
xticks([1:1:7])
grid on
ylabel('Frequency')
xlabel('Number of Synapses')
set(gca,'FontSize',18)
set(gca,'YGrid','off')

legend([p1,p2,p3,p4],{'Temporal+Hemilineage','Hemilineage','Temporal','Random'},'Location','best')

if savefigs == 1
    saveas(gcf,strcat(directory,'/','Network_Distances'),'svg')
    close all
else
end
%% Stats
[p_cohort_vs_hemi] = ranksum(tc_d,hemi_d)
[p_hemi_vs_temporal] = ranksum(temporal_d,hemi_d)
[p_temporal_vs_unrelated] = ranksum(temporal_d,ur_d)

%%
% find the number of one-hops and two-hops
one_hop_r = sum(ur_d==1)/numel(ur_d);
one_hop_h = sum(hemi_d==1)/numel(hemi_d);
one_hop_tc = sum(tc_d == 1)/numel(tc_d);
one_hop_t = sum(temporal_d == 1)/numel(temporal_d);

two_hops_r = sum(ur_d==2)/numel(ur_d);
two_hops_h = sum(hemi_d==2)/numel(hemi_d);
two_hops_tc = sum(tc_d ==2)/numel(tc_d);
two_hops_t = sum(temporal_d == 2)/numel(temporal_d);

three_hops_r = sum(ur_d>2)/numel(ur_d);
three_hops_h = sum(hemi_d>2)/numel(hemi_d);
three_hops_tc = sum(tc_d>2)/numel(tc_d);
three_hops_t = sum(temporal_d>2)/numel(temporal_d);

colors = [1,0,0;0,0,1;1,0,0;0,1,1]
figure;hold on
bg = barh([one_hop_r,one_hop_t,one_hop_h,one_hop_tc;two_hops_r,two_hops_t,two_hops_h,two_hops_tc;three_hops_r,three_hops_t,three_hops_h,three_hops_tc])
for k = 1:length(bg)
    bg(k).CData = repmat(colors(k,:),length(bg(k).CData))
end
yticks([1:3])
yticklabels({'Direct Connection','Two Synapses','> Two Synapses (max = 7)'})
xlabel('Fraction of Paths')
set(gca,'FontSize',18)
set(gca,'Ydir','reverse')
if savefigs == 1
    saveas(gcf,strcat(directory,'/','Network_Distances_barg'),'svg')
    close all
else
end
