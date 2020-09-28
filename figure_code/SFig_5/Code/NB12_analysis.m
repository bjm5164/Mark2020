%% First look at the synapse presynapse distribution of NB1-2 early born neurons compared to the rest of the lineage.
% Load dataset
load mark2020_neurons_temporal_cohorts.mat

% Get NB1-2 temporal cohorts
NB12 = nl(strcmp(an_in.Lineage(:),'NB1-2'));
NB12_an = an_in(strcmp(an_in.Lineage(:),'NB1-2'),:);
dorsal = NB12(NB12_an.DV_Index == 1);
dorsal_an = NB12_an(NB12_an.DV_Index == 1,:);

NB12_dorsal_early = nl(strcmp(an_in.Lineage(:),'NB1-2') & an_in.DV_Index == 1 & an_in.Temporal_Cohort == 1);
NB12_dorsal_mid = nl(strcmp(an_in.Lineage(:),'NB1-2') & an_in.DV_Index == 1 & an_in.Temporal_Cohort == 2);
NB12_dorsal_late1 = nl(strcmp(an_in.Lineage(:),'NB1-2') & an_in.DV_Index == 1 & an_in.Temporal_Cohort == 3);
NB12_dorsal_late2 = nl(strcmp(an_in.Lineage(:),'NB1-2') & an_in.DV_Index == 1 & an_in.Temporal_Cohort == 4);

% Plot synapse distributions
comb = [NB12_dorsal_early,NB12_dorsal_mid,NB12_dorsal_late1,NB12_dorsal_late2];
map = cbrewer('seq','Blues',12);
map([13 14],:) = [1 0 0;1,0,0];
map([15 16],:) = [.75 0 0; .75 0 0];
map([17 18 19 20],:) = [.5 0 0 ; .5 0 0 ; .5 0 0 ; .5 0 0 ];
plot_synapse_distributions(comb,map,1,.25)

%% Plot the neuropil distance of NB1-2 neurons
figure('pos',[1,1,200,1000]); hold on

for i = 1:length(comb)
    scatter(1,comb(i).Mean_Neuropil_Distance/1000,500,'o','MarkerEdgeColor','k','MarkerFaceColor',map(i,:));
end
xlim([0 2]);
xticks([]);
ylim([0 30]);
set (gca,'Ydir','reverse');
ylabel('Distance from Neuropil (µm)');
set(gca,'FontSize',24);
set(gcf,'Color','w');


%% Look at the average temporal cohort size compared to NB1-2 early

% Remove NB1-2 from the full neuron list
NB12D_index = boolean(strcmp(an_in.Lineage,'NB1-2') & an_in.DV_Index(:) == 1);
nl(NB12D_index) = [];
an_in(NB12D_index) = []; clear NB12D_index

% Get unique temporal cohort indices
tc_unique = unique(an_in(:,[2,4,5,6]),'rows');
[~,tc_index] = ismember(an_in(:,[2,4,5,6]),tc_unique,'rows')

% Compile counts
[C,ia,ic] = unique(tc_index);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];

% Boxplot of cohort sizes
figure; hold on
boxplot(a_counts)
scatter(1,6,'r','filled')
ylim([0 7])
xticklabels('')
ylabel('Cohort Size')
set(gca,'FontSize',18)

%% Look at similarity compared to other temporal cohorts
NB12_dorsal_early_presim = synapse_similarity_v2(NB12_dorsal_early,2000,3,[],1)
NB12_dorsal_early_postsim = synapse_similarity_v2(NB12_dorsal_early,2000,3,[],1)

presynaptic_similarity = NB12_dorsal_early_presim(boolean(triu(ones(size(NB12_dorsal_early_presim)),1)))
postsynaptic_similarity = NB12_dorsal_early_postsim(boolean(triu(ones(size(NB12_dorsal_early_postsim)),1)))

%% Compile data
NB1_2_analysis.Cohorts = {NB12_dorsal_early,NB12_dorsal_mid,NB12_dorsal_late1,NB12_dorsal_late2}
NB1_2_analysis.Synapse_Similarity = table(presynaptic_similarity,postsynaptic_similarity,'VariableNames',{'Presynapse_Similarity','Postsynapse_Similarity'})
NB1_2_analysis.Cohort_Sizes = a_counts

save(strcat(date(),'_','NB12_analysis'),'NB1_2_analysis')
