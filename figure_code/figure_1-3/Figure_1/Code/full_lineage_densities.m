load mark2020_neurons.mat
%% Parse lineages 
lineages = arrayfun(@(x) nl(an_in.Lineage_Index == x), 1:max(an_in.Lineage_Index),'UniformOutput',false)
lineages(cellfun(@isempty, lineages)) = []
lineage_legend = unique(an_in.Lineage)

%%
% Load neuropil bounds and rotate to correct for offset.
map = hsv(7)
load Neuropil_Mesh_Object.mat
xrot = 0;
yrot = -1;
zrot = -12;

NPM.vertices = rotate_points(NPM.vertices,zrot,3);
NPM.vertices = rotate_points(NPM.vertices,xrot,1);
NPM.vertices = rotate_points(NPM.vertices,yrot,2);


savefigs = input('Save Figures? 1:Yes 0:no')
if savefigs == 1
    directory = uigetdir
else
    directory = 0
end

% Plot 90% synapse density contours together.

figure('pos',[1,1,1200,1200],'rend','painters'); hold on
surfaces2D([NPM.vertices(NPM.vertices(:,3)> 105550 & NPM.vertices(:,3)< 144000,1),NPM.vertices(NPM.vertices(:,3)> 105550 & NPM.vertices(:,3)< 144000,2)],'k',.05,'-')    
for i = 1:length(lineages)
    synapse_density_contours(lineages{i},map(i,:),1,1,1,90) % Generate presynapse density contours for each hemilineage. Threshold is 75%, 1contour line. 
    axis off
    axis equal
    title(strcat('Presynaptic Density',lineage_legend{i}))
    if savefigs == 1
        saveas(gcf,string(strcat(directory,'/','All_Presynaptic','_Density.svg')),'svg')
    close all
    else
    end
end

% Plot individual 75% synapse density contours.
figure('pos',[1,1,1800,400],'rend','painters'); hold on
for i = 1:length(lineages)
    subplot(1,7,i)
    synapse_density_contours(lineages{i},'r',1,1,1,75) % Generate presynapse density contours for each hemilineage. Threshold is 75%, 1contour line. 
    surfaces2D([NPM.vertices(NPM.vertices(:,3)> 105550 & NPM.vertices(:,3)< 144000,1),NPM.vertices(NPM.vertices(:,3)> 105550 & NPM.vertices(:,3)< 144000,2)],'k',.05,'-')    
    axis off
    axis equal
    title(strcat('Presynaptic Density',lineage_legend{i}))
    
end
if savefigs == 1
        saveas(gcf,string(strcat(directory,'/','Individual_Presynaptic_Density.svg')),'svg')
    close all
    else
    end