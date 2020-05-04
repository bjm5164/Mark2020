%% Load neurons
load Hb
load Cas
load mark2020_neurons.mat
%% Get mean neurite length for each pair of neurons
for i = 1:2:length(Hb)
    hb_pairs = [Hb(i),Hb(i+1)];
    npd_l([i,i+1]) = hb_pairs(1).skeleton_data.Distance_To_Neuropil;
    npd_r([i,i+1]) = hb_pairs(2).skeleton_data.Distance_To_Neuropil;
    mean_npd([i,i+1]) = mean([npd_l(i),npd_r(i)]);
end

for i = 1:length(Hb)
    Hb(i).Mean_Neuropil_Distance = mean_npd(i);
end
clear npd_l and npd_r and mean_npd

for i = 1:2:length(Cas)
    cas_pairs = [Cas(i),Cas(i+1)];
    npd_l([i,i+1])= cas_pairs(1).skeleton_data.Distance_To_Neuropil;
    npd_r([i,i+1])= cas_pairs(2).skeleton_data.Distance_To_Neuropil;
    mean_npd([i,i+1]) = mean([npd_l(i),npd_r(i)]);
end


for i = 1:length(Cas)
    Cas(i).Mean_Neuropil_Distance = mean_npd(i);
end
clear npd_l and npd_r and mean_npd

%% Plot histogram of neurite lengths for known birthdates
edges = 0:2000:3e4

hb_neuropil_distances = [Hb(1:2:end).Mean_Neuropil_Distance];
cas_neuropil_distances = [Cas(1:2:end).Mean_Neuropil_Distance];

[N_hb,edges_cas] = histcounts(hb_neuropil_distances,edges,'Normalization','probability')
[N_cas,edges_cas] = histcounts(cas_neuropil_distances,edges,'Normalization','probability')
figure('pos',[0 0 815 1135]); hold on
%plot(edges(1:end-1),N_hb)
area(edges(1:end-1)/1000,N_hb,'FaceColor','c','FaceAlpha',.5,'EdgeColor','c','LineWidth',3)
%plot(edges(1:end-1),N_cas)
area(edges(1:end-1)/1000,N_cas,'FaceColor','r','FaceAlpha',.5,'EdgeColor','r','LineWidth',3)
view([270, 90])
xlabel('Cortex Neurite Length (µM)')
ylabel('Number of Neurons')
legend({'Hb','Cas'},'Location','best')
ax = gca ; ax.YAxisLocation = 'Right'
%set (gca,'Ydir','reverse')
ax.XAxisLocation = 'Top'
set(gcf,'Color','w')
set(gca,'FontSize',18)
set (gca,'Xdir','reverse')

% Set cohort edge limits
cohort_edge_lims(1) = 0
cohort_edge_lims(2) = mean(hb_neuropil_distances)+std(hb_neuropil_distances)
cohort_edge_lims(3) = mean(cas_neuropil_distances)-std(cas_neuropil_distances)
cohort_edge_lims(4) = mean(cas_neuropil_distances)

%final edge is the longest cortex neurite in the dataset
all_cortex_neurites = arrayfun(@(x) nl(x).skeleton_data.Distance_To_Neuropil,1:length(nl))
cohort_edge_lims(5) = max(all_cortex_neurites);
%% Generate temporal cohorts

[Neuron_List] = temporal_cohort_assignments(Neuron_List,'cell_body',0,1,1,cohort_edge_lims)
an = parse_annotations(Neuron_List,1)
nl = Neuron_List(an.DV_Index<2);
an_in = an(an.DV_Index<2,:);
an_in.Temporal_Cohort = transpose([nl(:).Temporal_Cohort])

%% Plot neurite distances colored by temporal cohort

lineages = arrayfun(@(x) nl(an_in.Lineage_Index == x), min(an_in.Lineage_Index):max(an_in.Lineage_Index),'UniformOutput',false);
lindex = arrayfun(@(x) an_in.Lineage(an_in.Lineage_Index == x), min(an_in.Lineage_Index):max(an_in.Lineage_Index),'UniformOutput',false);
