%% Organize the neuron list and calculate the birthdates.

% Generate the temporal cohorts.  This function has some options, but we just used the normalized cortex neurite length. 
Neuron_List = temporal_cohort_assignments(Neuron_List,'cell_body',1,1,1)

an = parse_annotations(Neuron_List,1); % Parse annotations
for i = 1:length(Neuron_List) % Rename neurons by paired names
    Neuron_List(i).Names_org = Neuron_List(i).Names;
    Neuron_List(i).Names = an.Names(i);
end
[~,si] = sort([Neuron_List(:).Names]); 
Neuron_List = Neuron_List(si); % Sort by paired names
an = an(si,:); clear si % Sort annotations by paired names
an.Temporal_Cohort = transpose([Neuron_List(:).Temporal_Cohort]);

nl = Neuron_List(an.DV_Index < 2); % Get just the interneurons
an_in = an(an.DV_Index<2,:); % Get just the interneuron annotations
save mark2020_neurons_organized



