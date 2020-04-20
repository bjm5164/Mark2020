function [Neuron_List] = temporal_cohort_assignments(Neuron_List,method,normalization,cutoff,bilateral)
% Generate the temporal cohort assignments for a neuron list.  There are
% two ways to do it:  
%               Cell_body: Distance clustering of the cell 
%               bodies of each lineage, and then assigns a temporal cohort 
%               based on the mean cortex neurite length for that cluster 
%               of cell bodies.
%               Cortex_neurite: Uses a kernel density estimate of cortex
%               neurite lenghts for a given lineage and looks for minima in
%               the density estimate to generate clusters of cells which
%               are then assigned a temporal cohort based on mean cortex
%               neurite length.
% Inputs:
%               Neuron_List: standard neuron list
%               method: either 'cell_body' or 'cortex_neurite' see above
%               for details.
%               
%               normalization:  1 for yes, 0 for no.  Normalization will
%               call normalized neurite distance to determine the cortex
%               neurite lenghts.  It essentially measures the with of the
%               cortex where the lineage is and normalizes the cortex
%               neurite length to that length.
%               
%               bandwidth: A distance in nm. This will change depending on
%               which method you use.  For cell_body, 5000 works well,
%               while for cortex_neurite, 1000-2000 works well.  Using a
%               very small bandwidth will cluster every cell into a unique
%               cluster which is equivalent to just using it's cortex
%               neurite length.
%               
%               Bilateral:  If the dataset is bilateral, it will average
%               the cortex neurite lengths for the l/r homologs.
% Outputs:      Neuron_List with temporal cohorts added.
%               

%% Parse lineages and hemilineages
an = parse_annotations(Neuron_List,1); % Parse annotations
for i = 1:length(Neuron_List) % Rename neurons by paired names
    Neuron_List(i).Names_org = Neuron_List(i).Names;
    Neuron_List(i).Names = an.Names(i);
end
[~,si] = sort([Neuron_List(:).Names]); 
Neuron_List = Neuron_List(si); % Sort by paired names
an = an(si,:);clear si % Sort annotations by paired names
in_index = find(an.DV_Index<2);
nl = Neuron_List(in_index); % Get just the interneurons
an_in = an(in_index,:); % Get just the interneuron annotations

lineage_indices = unique(an_in.Lineage_Index);
lineages = arrayfun(@(x) nl(an_in.Lineage_Index == x), min(an_in.Lineage_Index):max(an_in.Lineage_Index),'UniformOutput',false);
an_l = arrayfun(@(x) an_in(an_in.Lineage_Index == x,:), min(an_in.Lineage_Index):max(an_in.Lineage_Index),'UniformOutput',false);

% lineage_indices = unique([an_in.Lineage_Index,an_in.DV_Index],'rows')
% lineages = arrayfun(@(x) nl(an_in.Lineage_Index == lineage_indices(x,1) & an_in.DV_Index == lineage_indices(x,2)), 1:length(lineage_indices),'UniformOutput',false)
% an_l = arrayfun(@(x) an_in(an_in.Lineage_Index == lineage_indices(x,1) & an_in.DV_Index == lineage_indices(x,2),:), 1:length(lineage_indices),'UniformOutput',false)

%%

for i = 1:length(lineages)
    if bilateral == 1
        index_l = find(an_l{i}.Side_Index == 0); % index for left side neurons
        index_r = find(an_l{i}.Side_Index == 1); % index for right side neurons
        if normalization == 1
            npd_l = normalized_neurite_distance(lineages{i}(index_l),1); %Get cortex neurite lengths for left
            npd_r = normalized_neurite_distance(lineages{i}(index_r),1); %Get cortex neurite lengths for right
        else
            npd_l = arrayfun(@(x) lineages{i}(index_l(x)).skeleton_data.Distance_To_Neuropil,1:length(index_l)); %Get cortex neurite lengths for left
            npd_r = arrayfun(@(x) lineages{i}(index_r(x)).skeleton_data.Distance_To_Neuropil,1:length(index_r)); %Get cortex neurite lengths for right
        end
            
        npd = .5*(npd_l+npd_r); % Get mean cortex neurite length
        npd_bi = zeros(length(npd)*2,1);
        npd_bi(index_l) = npd;
        npd_bi(index_r) = npd;
    else
         % If not bilateral data, just get the cortex neurite lenghts
        if normalization == 1
            npd_bi = normalized_neurite_distance(lineages{i},1);
        else
            npd_bi = arrayfun(@(x) lineages{i}(x).skeleton_data.Distance_To_Neuropil,1:length(lineages{i}));
        end
        
    end

    for ii = 1:length(lineages{i})
        lineages{i}(ii).Mean_Neuropil_Distance = npd_bi(ii);
    end
    
    if contains(method,'cell_body')
    % Use cell body clustering
        cohorts{i} = get_temporal_cohorts_by_cell_body(lineages{i},cutoff,bilateral,normalization); % Get temporal cohorts by clustering cell body locations
        mean_c = cohorts{i}; % Take the mean for left and right homologs
        tc2 = mean_c;
    elseif contains(method,'cortex_neurite')
     % Use cortex neurite length
        cohorts{i} = get_temporal_cohorts_by_NPD(lineages{i},cutoff,bilateral,normalization); % Get temporal cohorts using a neurite length distance kernel
        mean_c = cohorts{i};
        tc2 = mean_c;
    else
        display('Using default cell body clustering')
        cohorts{i} = get_temporal_cohorts(lineages{i},cutoff,1); % Get temporal cohorts by clustering cell body locations
        mean_c = cohorts{i}; % Take the mean for left and right homologs
        tc2 = mean_c;
    end
    
    for kk = 1:length(lineages{i})
        lineages{i}(kk).Temporal_Cohort = tc2(kk);
    end
    
   
end
for i = 1:length(Neuron_List)
    Neuron_List(i).Temporal_Cohort = 0;
    Neuron_List(i).Mean_Neuropil_Distance = 0;
end

nl = cat(2,lineages{:})
for i = 1:length(lineages)
    nl(an_in.Lineage_Index == lineage_indices(i)) = lineages{i}
    %nl(ismember([an_in.Lineage_Index,an_in.DV_Index],lineage_indices(i,:),'rows'))
end


Neuron_List(in_index) = nl;