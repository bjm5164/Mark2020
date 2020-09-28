load mark2020_neurons_temporal_cohorts.mat
an.Temporal_Cohort = zeros(length(an.DV_Index),1)
an.Temporal_Cohort(an.Lineage_Index > 1 & an.Lineage_Index < 9) = an_in.Temporal_Cohort


neuron_deg = get_adjacency(Neuron_List,1);
neuron_deg(neuron_deg<.01) = 0;
neuron_deg(neuron_deg>0) = 1;
neuron_adj_ud = neuron_deg+neuron_deg';
neuron_deg(find(sum(neuron_adj_ud) == 0),:) = [];
neuron_deg(:,find(sum(neuron_adj_ud) == 0)) = [];
Neuron_List(find(sum(neuron_adj_ud) == 0)) = [];
an(find(sum(neuron_adj_ud) == 0),:) = []
side_lookup = containers.Map([0,1],{'Left','Right'});
hemi_lookup = containers.Map([0:3],{'Ventral','Dorsal','Motor','Sensory'});
lin_lookup = containers.Map(unique(an.Lineage_Index),unique(an.Lineage));
hemi_legend = unique(an.Lineage);
lineages = unique(an.Lineage);


[~,hemilineage_index] = ismember([an.Side_Index,an.DV_Index,an.Lineage_Index],unique([an.Side_Index,an.DV_Index,an.Lineage_Index],'rows'),'rows');

%% Sum_Connectivity
comb_adj = zeros(max(hemilineage_index));
for i = 1:max(hemilineage_index);
    hli_i = find(hemilineage_index == i);
    hemi_legend{i} = strcat(unique(an.Lineage(find(hemilineage_index == i))));
    for j = 1:max(hemilineage_index);
        hli_j = find(hemilineage_index == j);
        comb_adj(i,j) = sum(reshape(neuron_deg(hli_i,hli_j),numel(neuron_deg(hli_i,hli_j)),1));
        clear hli_j
    end
    clear hli_i
end

% Index values for removing the sensory inuts and motor outputs need to be
% adjusted depending on the dataset.
comb_adj_sansSM = comb_adj;
comb_adj_sansSM([13,27],:) = [];
comb_adj_sansSM(:,[14,28]) = [];

output_names = hemi_legend;
output_names([13,27]) = [];
input_names = hemi_legend;
input_names([14,28]) = [];


map = cbrewer('seq','Blues',64);
imagesc(comb_adj_sansSM,[0 10]); xticks(1:length(input_names)); 
xticklabels([input_names{:}]); xtickangle(90); 
yticks(1:length(output_names));yticklabels([output_names{:}]); 
axis xy;
colormap(map);
set(gca,'FontSize',18);
c = colorbar;
c.Label.String = '? Number of Connections (> 1% of total input)';
clear comb_adj_sansSM and input_names and output_names and map and c

adj = get_adjacency(Neuron_List,0);
comb_adj_frac_i = zeros(max(hemilineage_index));
for i = 1:max(hemilineage_index);
    hli_i = find(hemilineage_index == i);
    hemi_legend{i} = strcat(unique(an.Lineage(find(hemilineage_index == i)))); %,hemi_lookup(unique(an.DV_Index(find(hemilineage_index == i)))),' ',side_lookup(unique(an.Side_Index(find(hemilineage_index == i)))));
    for j = 1:max(hemilineage_index);
        hli_j = find(hemilineage_index == j);
        hemi_j = Neuron_List(hli_j);
        comb_adj_frac_i(i,j) = sum(reshape(adj(hli_i,hli_j),numel(adj(hli_i,hli_j)),1))/sum(arrayfun(@(x) length(hemi_j(x).Inputs.conid), 1:length(hemi_j)));
       
    end
end


s_map = cbrewer('qual','Set1',length(hemi_legend)/2);
double_map = ones(length(hemi_legend),3);
double_map(1:2:length(double_map),:) = s_map;
double_map(2:2:length(double_map),:) = s_map;

double_idx = ones(length(hemi_legend),1);
double_idx(1:2:end) = 1:1:length(hemi_legend)/2;
double_idx(2:2:end) = length(hemi_legend)/2+1:1:length(hemi_legend);

figure; hold on
re_ordered_adj_i = comb_adj_frac_i(double_idx,double_idx);
for i = 1:length(hemi_legend)
    bar(i,sum(re_ordered_adj_i(:,i)),'FaceColor',double_map(i,:));
end
xticks([1.5:2:length(hemi_legend)])
xticklabels([hemi_legend{double_idx(1:2:end)}])
ylabel('Fraction of Inputs')
set(gcf,'Color','w')
set(gca,'FontSize',18)

comb_adj_frac_o = zeros(max(hemilineage_index));
for i = 1:max(hemilineage_index);
    hli_i = find(hemilineage_index == i);
    hemi_legend{i} = strcat(unique(an.Lineage(find(hemilineage_index == i)))); 
    hemi_i = Neuron_List(hli_i);
    for j = 1:max(hemilineage_index);
        hli_j = find(hemilineage_index == j);
        hemi_j = Neuron_List(hli_j);
        comb_adj_frac_o(i,j) = sum(reshape(adj(hli_i,hli_j),numel(adj(hli_i,hli_j)),1))/sum(arrayfun(@(x) sum(hemi_i(x).Outputs.polyadics),1:length(hemi_i)));
       
    end
end


figure; hold on
re_ordered_adj_o = comb_adj_frac_o(double_idx,double_idx);
for i = 1:length(hemi_legend)
    bar(i,sum(re_ordered_adj_o(i,:)),'FaceColor',double_map(i,:))
end
xticks([1.5:2:length(hemi_legend)])
xticklabels([hemi_legend{double_idx(1:2:end)}])
ylabel('Fraction of Outputs')
set(gcf,'Color','w')
set(gca,'FontSize',18)
clear hemilineage_index

G = digraph(neuron_deg);
figure('Position',[1100 0 2000 2000],'rend','Painters'); hold on

pg = plot(G,'Layout','force','NodeColor',[.8,.8,.8],'EdgeAlpha',.5,'EdgeColor',[.8,.8,.8]);
mn_i_l = find(an.DV_Index == 2 & an.Side_Index == 0);
sn_i_l = find(an.DV_Index == 3 & an.Side_Index == 0);

sn_i_r = find(an.DV_Index == 3 & an.Side_Index == 1);
mn_i_r = find(an.DV_Index == 2 & an.Side_Index == 1);

for ii = 1:length(mn_i_l)
    highlight(pg,mn_i_l(ii),'NodeColor',[.5,.1,.8],'MarkerSize',10);
    highlight(pg,predecessors(G,mn_i_l(ii)),mn_i_l(ii),'EdgeColor',[.5,.1,.8]);
end
for ii = 1:length(sn_i_l)
    highlight(pg,sn_i_l(ii),'NodeColor','g','MarkerSize',10);
    highlight(pg,sn_i_l(ii), successors(G,sn_i_l(ii)), 'EdgeColor', 'g');
end

for ii = 1:length(mn_i_r)
    highlight(pg,mn_i_r(ii),'NodeColor',[.5,.1,.8]*.5,'MarkerSize',10);
    highlight(pg,predecessors(G,mn_i_r(ii)),mn_i_r(ii),'EdgeColor',[.5,.1,.8]*.5);
end
for ii = 1:length(sn_i_r)
    highlight(pg,sn_i_r(ii),'NodeColor',[0,.5,0],'MarkerSize',10);
    highlight(pg,sn_i_r(ii), successors(G,sn_i_r(ii)), 'EdgeColor', [0,.5,0]);
end


axis off; axis equal;
xlim([min(pg.XData),max(pg.XData)]);
ylim([min(pg.YData),max(pg.YData)]);
clear d and dist_edge;
set(gca,'FontSize',24);
pg.LineWidth = G.Edges.Weight;
axis off; 
set(gcf,'Color','w');


savefigs = input('Save Figures? 1:Yes 0:no')

if savefigs == 1
    directory = uigetdir
else
    directory = 0
end

map = hsv(length(unique(an.Lineage)));
G = digraph(neuron_deg);
for i = 1:max(an.Lineage_Index)
    n = find(an.Lineage_Index == i);
    figure('Position',[1100 0 2000 2000],'rend','Painters'); hold on
    pg = plot(G,'Layout','force','NodeColor',[.8,.8,.8],'EdgeAlpha',.5,'EdgeColor',[.8,.8,.8]);
    for ii = 1:length(n)     
        
        highlight(pg,n(ii),'NodeColor',map(i,:),'MarkerSize',20,'ArrowSize',10,'ArrowPosition',1);
        
    end
    xlim([min(pg.XData),max(pg.XData)])
    ylim([min(pg.YData),max(pg.YData)])
    clear d and dist_edge;
    title(lineages{i});
    set(gca,'FontSize',36)
    pg.LineWidth = G.Edges.Weight;
    axis off; axis equal;
    set(gcf,'Color','w');
    if savefigs == 1
        saveas(gcf,strcat(directory,'/',hemi_legend{i},'_network_graph'),'tiff')
        close all
    else
        %close all
     end
end

G = digraph(neuron_deg);
figure('Position',[1100 0 2000 2000],'rend','Painters'); 

for i = 1:max(an.Lineage_Index)
    subplot(3,3,i)
    n = find(an.Lineage_Index == i);
    n_hem = an.DV_Index(n);
    pg = plot(G,'Layout','force','NodeColor',[.9,.9,.9],'EdgeAlpha',.5,'EdgeColor',[.9,.9,.9],'ArrowPosition',1);
    for ii = 1:length(n)     
        if n_hem(ii) == 1 | n_hem(ii) == 2
            highlight(pg,n(ii),'NodeColor',[.5,.1,.8],'MarkerSize',8);
            highlight(pg,n(ii),successors(G,n(ii)),'EdgeColor',[.5,.1,.8],'LineWidth',1.5)
        elseif n_hem(ii) == 0 | n_hem(ii) == 3
            highlight(pg,n(ii),'NodeColor','g','MarkerSize',8);
            highlight(pg,n(ii),successors(G,n(ii)),'EdgeColor','g','LineWidth',1.5)
        
        else
        end
        
    end
    hold on
    if length(unique(n_hem))>1
        ax(2) = scatter(NaN,NaN,10,'g','MarkerFaceColor','g');
        ax(1) = scatter(NaN,NaN,10,[.5,.1,.8],'MarkerFaceColor',[.5,.1,.8]); 
        legend(ax,{'Dorsal','Ventral'},'Location','southeast');   
    elseif unique(n_hem) == 2
        ax(1) = scatter(NaN,NaN,10,[.5,.1,.8],'MarkerFaceColor',[.5,.1,.8]);
        legend(ax,{'Motor'},'Location','southeast');  
    elseif unique(n_hem) == 0
        ax(1) = scatter(NaN,NaN,10,'g','MarkerFaceColor','g');
        legend(ax,{'Ventral'},'Location','southeast');  
    else
        ax(1) = scatter(NaN,NaN,10,'g','MarkerFaceColor','g');
        legend(ax,{'Sensory'},'Location','southeast');  
    end
    clear ax
    clear d and dist_edge
    title(lineages{i})
    axis off; 
    axis equal
    xlim([min(pg.XData),max(pg.XData)]);
    ylim([min(pg.YData),max(pg.YData)]);
    set(gca,'FontSize',18)
    set(gcf,'Color','w')
    if savefigs == 1
        saveas(gcf,strcat(directory,'/',hemi_legend{i},'_network_graph_hl'),'tiff')
        close all
    else
        
     end

end


[~,hemilineage_index] = ismember([an.DV_Index,an.Lineage_Index],unique([an.DV_Index,an.Lineage_Index],'rows'),'rows');
temporal_color = [0,1,1; 0,0,1; 1,.65,.25 ;1,0,0]
for i = 1:max(hemilineage_index)-2
    n = find(hemilineage_index == i);
    distance = [0:7000:28000];
%     d = arrayfun(@(x) Neuron_List(n(x)).skeleton_data.Distance_To_Neuropil, 1:length(n));
%     d_mean = zeros(length(d),1);
%     d_mean(1:2:end) = .5*(d(1:2:end) + d(2:2:end));
%     d_mean(2:2:end) = d_mean(1:2:end);
%     d_total{i} = transpose([d(1:2:end);d(2:2:end)]);
%     dist_edge = discretize(d_mean,distance);
    dist_edge = [Neuron_List(n).Temporal_Cohort]
    %subplot(3,4,i)
    figure('Position',[1100 0 2000 2000],'rend','Painters'); hold on;
    pg = plot(G,'Layout','force','NodeColor',[.8 .8 .8],'EdgeAlpha',.5,'EdgeColor',[.8 .8 .8],'ArrowPosition',1);
    for ii = 1:length(n)
        highlight(pg,n(ii),'NodeColor',temporal_color(dist_edge(ii),:),'MarkerSize',8);
        if isempty(successors(G,n(ii)))
        else
        highlight(pg,n(ii),successors(G,n(ii)),'EdgeColor',temporal_color(dist_edge(ii),:),'LineWidth',1.5);
        end
    end
    axis equal
    xlim([min(pg.XData),max(pg.XData)]);
    ylim([min(pg.YData),max(pg.YData)]);
    %clear d and dist_edge;
    title(unique(an.Lineage(hemilineage_index == i)));
    set(gca,'FontSize',24);
    axis off; 
    set(gcf,'Color','w');
    
    
    
end

%% Measure network distances
 clearvars -except Neuron_List and an and neuron_deg 
 Gud = graph(neuron_deg+neuron_deg')
 Gud.Edges.Weight = ones(length(Gud.Edges.Weight),1)


ins = Neuron_List(an.DV_Index<2)

d_edge_total = zeros(length(Neuron_List),1);
d_edge_total(find(an.DV_Index<2)) = [ins(:).Temporal_Cohort]
for i = 1:5000
         %pick a lineage, and make sure it is a lineage with at least
         %two neurons in it. 
         w = 0;
         while w<1
                  
             n = randi([2,8]); % Pick a lineage
             s = randi([0,1]); % Pick a side
             l_index = find(an.Lineage_Index == n & an.Side_Index == s); % Get lineage index
             
             h = randi([0,1]); % Pick a hemilineage (Dorsal or ventral)
             h_index = find(an.Lineage_Index == n & an.DV_Index == h & an.Side_Index == s); % Get a hemilineage index
             
             if length(l_index) <2 | length(h_index) < 3  % Make sure the index at least contains two neurons
                 w = 0;
             else
                 % Pick a temporal cohort with at least two neurons in it
                 neurons = Neuron_List(h_index);
                 dist_edge = d_edge_total(h_index);
                 temptc_index = h_index(dist_edge == randi(max(dist_edge),1));
                 if length(temptc_index) < 2
                     w = 0;
                 else
                     w = 1;
                 end
              
             end
         end
         
         p_index_l = randperm(length(l_index),2);
         %test_l{i} = Neuron_List(l_index(p_index_l));
         net_dist_l(i) = distances(Gud,l_index(p_index_l(1)),l_index(p_index_l(2)));
         
         p_index_h = randperm(length(h_index),2);
         test_h{i} = Neuron_List(h_index(p_index_h)).SkIDs;
         net_dist_h(i) = distances(Gud,h_index(p_index_h(1)),h_index(p_index_h(2)));
         
         p_index_tc = randperm(length(temptc_index),2);
         test_tc{i} = [Neuron_List(temptc_index(p_index_tc)).SkIDs];
         net_dist_tc(i) = distances(Gud,temptc_index(p_index_tc(1)),temptc_index(p_index_tc(2)));
         
         % do the same for a temporal cohort across lineages
    
         t = randi([1,max(d_edge_total)]); % Pick a temporal cohort across lineages
         temporal_index = find(d_edge_total == t & an.Side_Index == s);
         
         p_index_t = randperm(length(temporal_index),2);
         %test_t{i} = Neuron_List(temporal_index(p_index_t));
         net_dist_t(i) = distances(Gud,temporal_index(p_index_t(1)),temporal_index(p_index_t(2)));
         %do the same for random interneurons
    
         in_index = find(an.DV_Index<2 & an.Side_Index == s ); % pick from interneurons from the same side
         p_index_r = randperm(length(in_index),2);
         %testr{i} = Neuron_List(in_index(p_index_r));
         net_dist_r(i) = distances(Gud,in_index(p_index_r(1)),in_index(p_index_r(2)));
         
        
    end

    
% Get rid of the unconnected measurements
net_dist_h(net_dist_h == inf) = []   
net_dist_l(net_dist_l == inf) = []   
net_dist_r(net_dist_r == inf) = []   
net_dist_t(net_dist_t == inf) = []  
net_dist_tc(net_dist_tc == inf) = []   

% Plot the synapse distances 
figure; hold on
for i = 1:7
    hop_r(i) = sum(net_dist_r==i)/numel(net_dist_r);
    hop_h(i) = sum(net_dist_h==i)/numel(net_dist_h);
    hop_tc(i) = sum(net_dist_tc == i)/numel(net_dist_tc);
    hop_t(i) = sum(net_dist_t == i)/numel(net_dist_t); 
    hop_l(i) = sum(net_dist_l == i)/numel(net_dist_l);
end
plot([1:7],hop_r,'LineWidth',2,'Color','k')
plot([1:7],hop_h,'LineWidth',2,'Color','r')
plot([1:7],hop_tc,'LineWidth',2,'Color','c')
plot([1:7],hop_t,'LineWidth',2,'Color','m')
plot([1:7],hop_l,'LineWidth',2,'Color','b')
xticks([1:1:7])
grid on
ylabel('Frequency')
xlabel('Number of Synapses')
set(gca,'FontSize',18)
set(gca,'YGrid','off')

% Plot synapses distances as a cumulative distribution
his_r = histcounts(net_dist_r,[.5:1:6.5],'Normalization','probability')
his_h = histcounts(net_dist_h,[.5:1:6.5],'Normalization','probability')
his_tc = histcounts(net_dist_tc,[.5:1:6.5],'Normalization','probability')
his_t = histcounts(net_dist_t,[.5:1:6.5],'Normalization','probability')

figure; hold on
p1 = plot([1:6],cumsum(his_tc),'c','LineWidth',5)
p2 = plot([1:6],cumsum(his_h),'r','LineWidth',5)
p3 = plot([1:6],cumsum(his_t),'m','LineWidth',5)
p4 = plot([1:6],cumsum(his_r),'k','LineWidth',5)


scatter([1:6],cumsum(his_r),1000,'.','k')
scatter([1:6],cumsum(his_h),1000,'.','r')
scatter([1:6],cumsum(his_tc),1000,'.','c')
scatter([1:6],cumsum(his_t),1000,'.','m')
xticks([1:1:7])
grid on
ylabel('Frequency')
xlabel('Number of Synapses')
set(gca,'FontSize',18)
set(gca,'YGrid','off')

legend([p1,p2,p3,p4],{'Temporal+Hemilineage','Hemilineage','Temporal','Random'},'Location','best')

% find the number of one-hops and two-hops
one_hop_r = sum(net_dist_r==1)/numel(net_dist_r);
one_hop_h = sum(net_dist_h==1)/numel(net_dist_h);
one_hop_tc = sum(net_dist_tc == 1)/numel(net_dist_tc);
one_hop_t = sum(net_dist_t == 1)/numel(net_dist_t);
one_hop_l = sum(net_dist_l == 1)/numel(net_dist_l);

two_hops_r = sum(net_dist_r==2)/numel(net_dist_r);
two_hops_h = sum(net_dist_h==2)/numel(net_dist_h);
two_hops_tc = sum(net_dist_tc ==2)/numel(net_dist_tc);
two_hops_t = sum(net_dist_t == 2)/numel(net_dist_t);
two_hops_l = sum(net_dist_l == 2)/numel(net_dist_l);

three_hops_r = sum(net_dist_r>2)/numel(net_dist_r);
three_hops_h = sum(net_dist_h>2)/numel(net_dist_h);
three_hops_tc = sum(net_dist_tc>2)/numel(net_dist_tc);
three_hops_t = sum(net_dist_t>2)/numel(net_dist_t);
three_hops_l = sum(net_dist_l>2)/numel(net_dist_l);

colors = [1,0,0;0,0,1;1,0,0;0,1,1]
figure;hold on
bg = barh([one_hop_r,one_hop_t,one_hop_h,one_hop_tc;two_hops_r,two_hops_t,two_hops_h,two_hops_tc;three_hops_r,three_hops_t,three_hops_h,three_hops_tc])
for k = 1:length(bg)
    bg(k).CData = repmat(colors(k,:),length(bg(k).CData))
end
yticks([1:3])
yticklabels({'Direct Connection','Two Synapses','> Two Synapses (max = 7)'})
xlabel('Fraction of Paths')
set(gca,'FontSize',18)
set(gca,'Ydir','reverse')
[h,p] = ranksum(net_dist_h,net_dist_tc)






