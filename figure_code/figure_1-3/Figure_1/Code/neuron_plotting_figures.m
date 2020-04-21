% This will generate the plots for all interneurons colored by either
% hemilineage or lineage.  

load mark2020_neurons.mat
an = parse_annotations(Neuron_List,1)
%% Plot all neurons colored by lineage
index = [2:8]
map = hsv(length(index))
figure; hold on
for i = 1:length(index)
    plot_neurons(Neuron_List(an.Lineage_Index == index(i)),map(i,:),1,3,1,0)
end
view([190 90]);
camlight ;
axis off;
axis equal;
set(gcf,'Color','w');


%% Plot Hemilineages
plot_neurons(Neuron_List(an.DV_Index == 1),'r',1,3,1,0)
plot_neurons(Neuron_List(an.DV_Index == 0),'b',1,3,1,0)
view([190 90]);
camlight ;
axis off;
axis equal;
set(gcf,'Color','w');