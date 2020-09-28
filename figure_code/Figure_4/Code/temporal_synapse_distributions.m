load mark2020_neurons_temporal_cohorts.mat

% Remove NB1-2D from this analysis.
NB12D_index = boolean(strcmp(an_in.Lineage,'NB1-2') & an_in.DV_Index(:) == 1);
nl(NB12D_index) = [];
an_in(NB12D_index,:) = []; clear NB12D_index

%Separate into dorsal and ventral hemilineages
dorsal_neurons = nl(an_in.DV_Index == 1)
an_dorsal = an_in(an_in.DV_Index == 1,:)
ventral_neurons = nl(an_in.DV_Index == 0)
an_ventral = an_in(an_in.DV_Index == 0,:)

%%
map =  [.25,.9,.9; .25,.25,.9; .9,.55,.25 ;.9,.25,.25]

for i = 1:4
    tc_d{i} = dorsal_neurons(an_dorsal.Temporal_Cohort == i)
    tc_v{i} = ventral_neurons(an_ventral.Temporal_Cohort == i)
end
    
plot_synapse_distributions(tc_d,map,1,.25)
plot_synapse_distributions(tc_d,map,2,.25)

plot_synapse_distributions(tc_v,map,1,.25)
plot_synapse_distributions(tc_v,map,2,.25)