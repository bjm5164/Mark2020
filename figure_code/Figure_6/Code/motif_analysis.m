load mark2020_neurons_temporal_cohorts.mat
an.Temporal_Cohort = zeros(length(an.DV_Index),1)
an.Temporal_Cohort(an.Lineage_Index > 1 & an.Lineage_Index < 9) = an_in.Temporal_Cohort


neuron_deg = get_adjacency(Neuron_List,1);
neuron_deg(neuron_deg<.01) = 0;
neuron_deg(neuron_deg>0) = 1;
neuron_adj_ud = neuron_deg+neuron_deg';
unconnected_index = sum(neuron_adj_ud) == 0;
neuron_deg(unconnected_index,:) = [];
neuron_deg(:,unconnected_index) = []; 
Neuron_List(unconnected_index) = [];
an(unconnected_index,:) = []; clear unconnected_index
side_lookup = containers.Map([0,1],{'Left','Right'});
hemi_lookup = containers.Map([0:3],{'Ventral','Dorsal','Motor','Sensory'});
lin_lookup = containers.Map(unique(an.Lineage_Index),unique(an.Lineage));
hemi_legend = unique(an.Lineage);
lineages = unique(an.Lineage);


[~,hemilineage_index] = ismember([an.Side_Index,an.DV_Index,an.Lineage_Index],unique([an.Side_Index,an.DV_Index,an.Lineage_Index],'rows'),'rows');

%%
recurrent = zeros(length(neuron_deg))
common_in = zeros(length(neuron_deg))
common_out = zeros(length(neuron_deg))
for i = 1:length(neuron_deg); 
    rec = (neuron_deg(:,i) + neuron_deg(i,:)');
    recurrent(i,rec==2) = 1;
    for j = 1:length(neuron_deg);
   
        common_out = neuron_deg(i,:) + neuron_deg(j,:);
        common_in = neuron_deg(:,i) + neuron_deg(:,j);
        
        
        shared_in(i,j) = sum(common_in == 2);
        shared_out(i,j) = sum(common_out == 2);
        
        clear common_in and common_out
        
    end
    
end
%%

% Get network distances for temporal cohort related neurons
temporal_cohort_options = unique(an(:,[2,4,6,5]),'rows')
temporal_cohort_options(temporal_cohort_options.DV_Index>1,:) = []
% Remove NB1-2 from temporal cohort analysis
temporal_cohort_options(ismember(table2array(temporal_cohort_options) , [1 2 1 1], 'rows')...
                        | ismember(table2array(temporal_cohort_options) , [1 2 1 0], 'rows'),:) = []

for i = 1:height(temporal_cohort_options)
    related_index = ismember(an(:,[2,4,6,5]),temporal_cohort_options(i,:));
    temp_distances = recurrent(related_index,related_index);
    
    %Take only the single comparisons
    temp_distances = temp_distances(boolean(triu(ones(size(temp_distances)),1)))
    
    recurrent_tc{i} = temp_distances; clear temp_distances
end
%%

shared_in(boolean(eye(size(shared_in)))) = 0
shared_out(boolean(eye(size(shared_out)))) = 0

%%
[ii,jj] = ind2sub(size(shared_out),find(shared_out > 1))
table(transpose([Neuron_List(ii).Names]),transpose([Neuron_List(jj).Names]))
table(table2array(an(ii,[2,4,5,6])),table2array(an(jj,[2,4,5,6])))

%%

temporal_cohort_options = unique(an(:,[2,4,6,5]),'rows')
temporal_cohort_options(temporal_cohort_options.DV_Index>1,:) = []

temporal_cohort_options(temporal_cohort_options.DV_Index==0,:) = []
% Remove NB1-2 from temporal cohort analysis
temporal_cohort_options(ismember(table2array(temporal_cohort_options) , [1 2 1 1], 'rows')...
                        | ismember(table2array(temporal_cohort_options) , [1 2 1 0], 'rows'),:) = []
                    
                    
for i = 1:height(temporal_cohort_options)
    index = ismember(an(:,[2,4,6,5]),temporal_cohort_options(i,:));
    if sum(index) > 1

        group_shared_in = shared_in(index,index)
        tc_shared_in{i} = group_shared_in(boolean(triu(ones(size(group_shared_in)),1)));
        
        group_shared_out = shared_out(index,index)
        tc_shared_out{i} = group_shared_out(boolean(triu(ones(size(group_shared_out)),1)));
    else
     
    end
end
%%
tc_in = cat(1,tc_shared_in{:})
tc_out = cat(1,tc_shared_out{:})

 
figure; hold on
histogram(tc_out_d,[0:1:5],'Normalization','probability','FaceColor','k')
histogram(tc_out,[0:1:5],'Normalization','probability','FaceColor','r')
