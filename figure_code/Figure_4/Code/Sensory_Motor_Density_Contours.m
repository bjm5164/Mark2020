load PMNs.mat
load post_sensory.mat
load m19_neurons_original.mat
MNs = Neuron_List(an.Lineage_Index == 1)
sensory = Neuron_List(an.Lineage_Index == 9)
%% Make synapse distributions for Sensory and Motor Neurons.  Synapses are post-synaptic for motor neurons and pre-synaptic for sensory neurons.
figure('pos',[1,1,1800,1200]); hold on
Synapse_Distributions(sensory,'g',1,.25,0)
Synapse_Distributions(MNs,[.5 0 .75],2,.25,1)

%% Make synapse distributions for Post-sensory and pre-motor neurons.  Synapses are post-synaptic for pre-motors and pre-synaptic for post-sensory neurons.  

figure('pos',[1,1,1800,1200]); hold on

Synapse_Distributions(post_sensory,'g',1,.25,0)
Synapse_Distributions(PMNs,[.5 0 .75],2,.15,0)

%% Make synapse density contours for Post-sensory and pre-motor neurons. Given the limitations of contour plots in MATLAB, these have to be edited elsewhere.  
%To make the publication figure, the colormaps were used to assign
%different alpha values such that lower density contours are more
%transparent.  Then all of the contours are changed to the same color such
%that density contour changes are represented by decreasing alpha values,
%not color changes.
load Neuropil_Mesh_Object.mat
xrot = 0;
yrot = -1;
zrot = -12;

NPM.v = rotate_pointsV2(NPM.v,zrot,3);
NPM.v = rotate_pointsV2(NPM.v,xrot,1);
NPM.v = rotate_pointsV2(NPM.v,yrot,2);


contour_num = 1
figure('pos',[1,1,1800,1200],'rend','painters'); hold on
surfaces([NPM.v(NPM.v(:,3)>.9e5,1),NPM.v(NPM.v(:,3)>.9e5,2)],'k',.05,'-')
synapse_density_contours_v2(PMNs,'b',3,contour_num,1,70)
synapse_density_contours_v2(post_sensory,'g',3,contour_num,1,70)
axis off
axis equal