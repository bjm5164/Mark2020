savefigs = input('Save Figures? 1:Yes 0:no')
if savefigs == 1
    directory = uigetdir
else
    directory = 0
end

%% Plot synapse distributions for each ventral hemilineage
ventral_index = unique([an.Lineage_Index(an.DV_Index == 0)])
map = hsv(max(ventral_index))
map(end-1,:) = [1 0 1]
map(end-2,:) = [0.3 0 1.0000]
map(end-3,:) = [0 .6 1]


ventral_legend = unique([an.Lineage(an.DV_Index == 0)])
ventral_hl = arrayfun(@(x) Neuron_List(an.Lineage_Index == ventral_index(x,1) & an.DV_Index == 0 & an.Side_Index == 0), 1:length(ventral_index),'UniformOutput',false)

plot_synapse_distributions(ventral_hl,map(ventral_index-1,:),1,.25)
ax = subplot(3,7,[2 3 4 9 10 11])
legend([ax.Children(1:end-1)],ventral_legend,'Location','south')
clear ax
if savefigs == 1
        saveas(gcf,strcat(directory,'/','Ventral_hl_dist_pre'),'svg')
        close all
    else
end


plot_synapse_distributions(ventral_hl,map(ventral_index-1,:),2,.5)
ax = subplot(3,7,[2 3 4 9 10 11])
legend([ax.Children(1:end-1)],ventral_legend,'Location','south')
clear ax
if savefigs == 1
        saveas(gcf,strcat(directory,'/','Ventral_hl_dist_post'),'svg')
        close all
    else
end


%% Plot synapse distributions for each dorsal hemilineage     
dorsal_index = unique([an.Lineage_Index(an.DV_Index == 1)])
dorsal_legend = unique([an.Lineage(an.DV_Index == 1)])
dorsal_hl = arrayfun(@(x) Neuron_List(an.Lineage_Index == dorsal_index(x,1) & an.DV_Index == 1 & an.Side_Index == 0), 1:length(dorsal_index),'UniformOutput',false)

plot_synapse_distributions(dorsal_hl,map(dorsal_index-1,:),1,.25)
ax = subplot(3,7,[2 3 4 9 10 11])
legend([ax.Children(1:end-1)],dorsal_legend,'Location','south')
clear ax
if savefigs == 1
        saveas(gcf,strcat(directory,'/','Dorsal_hl_dist_pre'),'svg')
        close all
    else
end

plot_synapse_distributions(dorsal_hl,map(dorsal_index-1,:),2,.5)
ax = subplot(3,7,[2 3 4 9 10 11])
legend([ax.Children(1:end-1)],dorsal_legend,'Location','south')
clear ax
if savefigs == 1
        saveas(gcf,strcat(directory,'/','Dorsal_hl_dist_post'),'svg')
        close all
    else
end

