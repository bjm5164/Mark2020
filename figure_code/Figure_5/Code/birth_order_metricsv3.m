function [output_vars] = birth_order_metricsv3(nl,direction,suppress_figs,savefigs,directory)
if direction == 1
    in_out_dist = 'presynaptic';
    in_out_dvs = 'presynaptic';
    in_out_clust = 'presynaptic';
    io_dist = 1;
    io_dvs = 1;
elseif direction == 2
    in_out_dist = 'postsynaptic';
    in_out_dvs = 'postsynaptic';
    in_out_clust = 'postsynaptic';
    io_dist = 2;
    io_dvs = 2;
elseif direction == 3
    in_out_dist = 'presynaptic';
    in_out_dvs = 'combinedsynaptic';
    in_out_clust = 'combinedsynaptic';
    io_dist = 1;
    io_dvs = 3;
else error('Incorrect Direction')
end
an_nl = parse_annotations(nl,1);

if unique(an_nl.DV_Index) == 1
    dv = 'Dorsal'
elseif unique(an_nl.DV_Index == 0)
    dv = 'Ventral'
else
    warning('Not a hemilineage')
end
lineage = unique(an_nl.Lineage)
distance = [0,.39,.64,.85,1.2];
map = [.25,.9,.9; .25,.25,.9; .9,.55,.25 ;.9,.25,.25]

edges = [nl(:).Temporal_Cohort]
d = [nl(:).Mean_Neuropil_Distance]

%% synapse_distributions

% f1 = figure('pos',[1975 153 1375 715],'rend','painters'); hold on
% set(0, 'CurrentFigure', f1)

for i = 1:max(edges)
    cohort_list{i} = nl(edges == i)
end
tmap = map(~cellfun('isempty',cohort_list),:)
cohort_list(cellfun('isempty',cohort_list)) = []

plot_synapse_distributions(cohort_list,tmap,io_dist,.5)
clear cohort_list and tmap
% Synapse_Distributions(nl(1),map(edges(1),:),io_dist,.5,1)
% for i = 2:length(nl)
%     set(0, 'CurrentFigure', f1)
%     Synapse_Distributions(nl(i),map(edges(i),:),io_dist,.5,0)
% end

if savefigs == 1
    saveas(gcf,strcat(directory,'/',char(lineage),dv,in_out_dist,'_Distributions'),'svg')
    close all
else
end


%% Skeletons
figure('pos',[1975 153 1375 1375]); hold on
for i = 1:length(nl)
    plot_neurons(nl(i),map(edges(i),:),1,3,1,0);
end

    
camlight;
axis off; 
axis equal;
set(gcf,'Color','w');
view([190 90]);
c = colorbar;
c.Label.String = 'Distance from Neuropil (nm)';
c.Ticks = distance;
caxis([0 30000]);
colormap(map);
set(gca,'FontSize',18)    
if savefigs == 1
        print(gcf,strcat(directory,'/',char(lineage),dv,in_out_dist,'_Skeletons'),'-dtiff','-r300')
    close all
    else
end

    
%% Dotplot Neuropil Distnace
figure('pos',[1,1,200,1000]); hold on

%errorbar(mean(d),std(d),'LineWidth',2,'Color','k')
plot([.9,1.1],[mean(d),mean(d)],'k','LineWidth',2)
for i = 1:length(nl)
    scatter(1,d(i),500,'o','MarkerEdgeColor','k','MarkerFaceColor',map(edges(i),:))
end
xlim([0 2])
xticks([])

set (gca,'Ydir','reverse')
ylabel('Distance from Neuropil (nm)')
set(gca,'FontSize',18)
set(gcf,'Color','w')

if savefigs == 1
    saveas(gcf,strcat(directory,'/',char(lineage),dv,'_Neuropil_Distances'),'svg')
    close all
else
end



%% Cluster synapses and color labels by birth order
d_half = d(1:2:end)
edges_half = edges(1:2:end)

clear empty_ind
for i = 1:length(nl)
    if isempty(nl(i).Inputs.treenodeID) == 1
        empty_ind(i) = 0
    elseif isempty(nl(i).Outputs.treenodeID) == 1
        empty_ind(i) = 0
    else
        empty_ind(i) = 1
    end
end
nl = nl(boolean(empty_ind))
edges = edges(boolean(empty_ind))

if direction < 3
    sim_mat_l = synapse_similarity_v2(nl(an_nl.Side_Index == 0),4000,3,[],direction)
    sim_mat_r = synapse_similarity_v2(nl(an_nl.Side_Index == 1),4000,3,[],direction)
    sim_mat = .5*(sim_mat_l + sim_mat_r)
elseif direction == 3
    presim_mat_l = synapse_similarity_v2(nl(an_nl.Side_Index == 0),4000,3,[],1)
    presim_mat_r = synapse_similarity_v2(nl(an_nl.Side_Index == 1),4000,3,1)
    presim_mat = .5*(presim_mat_l + presim_mat_r)
    
    postsim_mat_l = synapse_similarity_v2(nl(an_nl.Side_Index == 0),4000,3,[],2)
    postsim_mat_r = synapse_similarity_v2(nl(an_nl.Side_Index == 1),4000,3,[],2)
    postsim_mat = .5*(postsim_mat_l + postsim_mat_r)
    sim_mat = (presim_mat+postsim_mat)*.5;
else error('Wrong Direction')
end

[h1 t1 perm1 Z] = Synapse_Distance_Clustering_v2(sim_mat,[nl(an_nl.Side_Index == 0).Names],0,0)
reordered_edges = edges_half(perm1);
subplot(3,4,[1 5 9]);
ax = get(gca)
lab = ax.YAxis.TickLabels

for i = 1:length(lab)
    end_name_ind = strfind(lab{i},'_')
    lab{i}(end_name_ind:end) = []
end

loc = ax.YAxis.TickValues;
for k = 1:numel(lab) % for every type of lable
    x = repelem(ax.XAxis.Limits(1)-0.01,1); % make an x position vector
    text(x,loc(k),lab(k),'Color',map(reordered_edges(k),:)); % place this lable at the same locations with a distinct color:
    ax.YAxis.TickLabels = []; % remove the original labels
    ax.YAxis.TickLabels = repelem('  ',max(cellfun(@numel,lab))); % replace the original labels with white space, to keep the axes position:
    set(gca,'FontSize',18)
end
if savefigs == 1
    saveas(gcf,strcat(directory,'/',char(lineage),dv,in_out_clust,'_Similarity_Clustering'),'svg')
    close all
else
end
%% Distance from Neuropil vs Similarity to other neurons in the hemilineage
figure; hold on
for i = 1:length(sim_mat)
    scatter(d_half(i),mean(sim_mat(i,:)),200,'MarkerEdgeColor','k','MarkerFaceColor',map(edges_half(i),:),'MarkerEdgeAlpha',.5)
    dist_vs_sim(i,:) = [d_half(i),mean(sim_mat(i,:))]
end
if savefigs == 1
    saveas(gcf,strcat(directory,'/',char(lineage),dv,in_out_dvs,'_Dist_vs_Sim'),'svg')
    close all
else
end

% Measure the closest (youngest) neuron of the pair
for i = 1:length(d)
    for j = 1:length(d)
        min_distance_matrix(i,j) = min([d(i),d(j)])
    end
end
        

%% Calculate Distance Differences
left_ind = find(an_nl.Side_Index == 0)
right_ind = find(an_nl.Side_Index == 1)

for i = 1:length(left_ind)
    for j = 1:length(left_ind)
        birth_order_distance_l(i,j) = pdist2(nl(left_ind(i)).Soma_Coords,nl(left_ind(j)).Soma_Coords);
        birth_order_distance_r(i,j) = pdist2(nl(right_ind(i)).Soma_Coords,nl(right_ind(j)).Soma_Coords);
    end
end
mean_birth_order_distance = .5*(birth_order_distance_l+birth_order_distance_r);
birth_order_distance = 1-(mean_birth_order_distance./max(mean_birth_order_distance(:)))
%% Quantify similarity
clear lin_sim and lin_bo
lin_sim = sim_mat(:);
lin_bo = birth_order_distance(:);

si = find(lin_sim == 1);
lin_sim(si) = []
lin_bo(si) = []
clear si

figure('Position',[1100 0 2000 900],'rend','Painters'); hold on
scatter(lin_bo,lin_sim,200,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerFaceAlpha',.5);
%[f_pre,gof,output]  = fit(lin_bo,lin_sim,'gauss1')
f = fitlm(lin_bo,lin_sim)
p = plot(f)


xlabel('Birth Order Similarity Metric')
ylabel(strcat(in_out_dvs,'Similarity'))
set(gca,'FontSize',18)
set(gcf,'Color','w')
xlim([0,1])
if savefigs == 1
    saveas(gcf,strcat(directory,'/',char(lineage),dv,'_birth_order_vs',in_out_dvs,'_similarity'),'svg')
    close all
else
end

%% Quantify similarity 

for i = 1:100
    % get a pair of neurons from the same temporal cohort
    
    % Pick a temporal cohort, and make sure it has more than one neuron in
    % it
    j = 0
    while j<1
        t_d = edges_half(randi(length(edges_half),1));
        if sum(edges_half==t_d) >= 2
            j = 1;
        else
            j = 0;
        end
    end
    % Get the synapse similarity matrix for the temporal cohort
    temp_sim_d = sim_mat(edges_half==t_d,edges_half==t_d); clear t_d
    
    % Pick a pair from that temporal cohort and get their similarity, then
    % make sure the pair isn't self-self.
    jj = 0
    while jj<1  
        pair_td = randperm(length(temp_sim_d),2);
        if pair_td(1) == pair_td(2)
            jj = 0
        else
            jj = 1
        end
        sim_td(i) = temp_sim_d(pair_td(1),pair_td(2)); clear pair_td
    end
    %Pick two random neurons from the hemilineage, and get their similarity
    k = 0
    while k<1
        t_rand = randi(length(edges_half),1,2);
        if t_rand(1) == t_rand(2)
            k = 0;
        else
            k = 1;
        end
    end
    sim_rand(i) = sim_mat(t_rand(1),t_rand(2));    
end

temporal_similarity = [sim_td(:),sim_rand(:)];
% figure; hold on
% scatter(ones(100,1),sim_td,'jitter','on')
% scatter(ones(100,1)+1,sim_rand,'jitter','on')
% xlim([0,3])
if suppress_figs == 1
    close all
else
end

output_vars.dist_vs_sim = dist_vs_sim
output_vars.sim_mat = sim_mat
output_vars.birth_order_distance = birth_order_distance
output_vars.model = f
output_vars.temporal_similarity = temporal_similarity
output_vars.min_distance_matrix = min_distance_matrix

end


