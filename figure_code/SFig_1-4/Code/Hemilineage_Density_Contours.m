load mark2020_neurons.mat
% Parse lineages and hemilineages
hemilineage_index = unique([an.Lineage_Index,an.DV_Index],'rows')

lineages = arrayfun(@(x) Neuron_List(an.Lineage_Index == x), 1:max(an.Lineage_Index),'UniformOutput',false)
hemilineages = arrayfun(@(x) Neuron_List(an.Lineage_Index == hemilineage_index(x,1) & an.DV_Index == hemilineage_index(x,2)), 1:length(hemilineage_index),'UniformOutput',false)

for i = 1:length(hemilineage_index)
    if hemilineage_index(i,2) == 0
      hemi_legend{i} = strcat(unique(an.Lineage(an.Lineage_Index == hemilineage_index(i,1))),'Ventral')
    elseif hemilineage_index(i,2) == 1
      hemi_legend{i} = strcat(unique(an.Lineage(an.Lineage_Index == hemilineage_index(i,1))),'Dorsal')
    else
      hemi_legend{i} = unique(an.Lineage(an.Lineage_Index == hemilineage_index(i,1)))
    end
end

%%
% Load neuropil bounds and rotate to correct for offset.
load Neuropil_Mesh_Object.mat
xrot = 0;
yrot = -1;
zrot = -12;

NPM.vert = rotate_pointsV2(NPM.vert,zrot,3);
NPM.vert = rotate_pointsV2(NPM.vert,xrot,1);
NPM.vert = rotate_pointsV2(NPM.vert,yrot,2);

savefigs = input('Save Figures? 1:Yes 0:no')

if savefigs == 1
    directory = uigetdir
else
    directory = 0
end

for i = 2:length(hemilineages)-1
    figure('pos',[1,1,1800,1200],'rend','painters'); hold on
    synapse_density_contours_v2(hemilineages{i},'r',1,10,1,60) % Generate presynapse density contours for each hemilineage. Threshold is 60%, 10 contour lines. 
    surfaces([NPM.vert(NPM.v(:,3)> 105550 & NPM.vert(:,3)< 144000,1),NPM.vert(NPM.v(:,3)> 105550 & NPM.vert(:,3)< 144000,2)],'k',.05,'-')    
    axis off
    axis equal
    title(strcat('Presynaptic Density',hemi_legend{i}))
    if savefigs == 1
        saveas(gcf,string(strcat(directory,'/',hemi_legend{i},'_Presynaptic','_Density.svg')),'svg')
    close all
    else
    end
end

for i = 2:length(hemilineages)-1
    figure('pos',[1,1,1800,1200],'rend','painters'); hold on
    synapse_density_contours_v2(hemilineages{i},'r',2,10,1,60) % Generate postsynapse density contours for each hemilineage. Threshold is 60%, 10 contour lines. 
    surfaces([NPM.vert(NPM.vert(:,3)>.9e5,1),NPM.vert(NPM.vert(:,3)>.9e5,2)],'k',.05,'-')
    axis off
    axis equal
    title(strcat('Postsynaptic Density',hemi_legend{i}))
    if savefigs == 1
        saveas(gcf,string(strcat(directory,'/',hemi_legend{i},'_Postsynaptic','_Density.svg')),'svg')
    close all
    else
    end
end




