load mark2020_neurons_temporal_cohorts.mat
savefigs = input('Save Figures? 1:Yes 0:no')
if savefigs == 1
    directory = uigetdir
else
    directory = 0

    
end

nl(1:24) = []
an_in(1:24,:) = []
%% Parse lineages and hemilineages

hemilineage_index = unique([an_in.Lineage_Index,an_in.DV_Index],'rows')

hemilineages = arrayfun(@(x) nl(an_in.Lineage_Index == hemilineage_index(x,1) & an_in.DV_Index == hemilineage_index(x,2)), 1:length(hemilineage_index),'UniformOutput',false)

for i = 1:length(hemilineage_index)
    if hemilineage_index(i,2) == 0
      hemi_legend{i} = strcat(unique(an_in.Lineage(an_in.Lineage_Index == hemilineage_index(i,1))),'Ventral')
    elseif hemilineage_index(i,2) == 1
      hemi_legend{i} = strcat(unique(an_in.Lineage(an_in.Lineage_Index == hemilineage_index(i,1))),'Dorsal')
    else
      hemi_legend{i} = unique(an_in.Lineage(an_in.Lineage_Index == hemilineage_index(i,1)))
    end
end 


dorsal_hemilineages = hemilineages(hemilineage_index(:,2) == 1);
dorsal_hemilineages{1} = []
ventral_hemilineages = hemilineages(hemilineage_index(:,2) == 0);
clearvars -except dorsal_hemilineages and ventral_hemilineages and hemilineages and Neuron_List and savefigs and directory an an and hemi_legend and nl and an_in
%% Birth Order Analysis of Dorsal Hemilineage Presynapses


% Calculate metrics
for i = 1:length(dorsal_hemilineages)
    if length(dorsal_hemilineages{i}) <= 2
    else
    temp_stats = birth_order_metricsv3(dorsal_hemilineages{i},1,0,savefigs,directory)
    dorsal_pre(i).temporal_similarity = temp_stats.temporal_similarity;
    dorsal_pre(i).birth_order_sim = temp_stats.birth_order_distance;
    dorsal_pre(i).sim_mat = temp_stats.sim_mat
    dorsal_pre(i).dvs = temp_stats.dist_vs_sim
    dorsal_pre(i).model = temp_stats.model
    dorsal_pre(i).neurite_lengths = [temp_stats.min_distance_matrix(:),temp_stats.mean_distance_matrix(:)]
    clear temp_stats
    end
end

% Compile synapse similarity, birthorder similarity, and average age
clear dorsal_presynaptic_similarity and dorsal_birthorder_similarity and dorsal_neurite_lengths
dorsal_presynaptic_similarity = dorsal_pre(1).sim_mat(:)
dorsal_birthorder_similarity = dorsal_pre(1).birth_order_sim(:)
dorsal_neurite_lengths = dorsal_pre(1).neurite_lengths;
for i = 2:length(dorsal_pre)
    dorsal_presynaptic_similarity = [dorsal_presynaptic_similarity; dorsal_pre(i).sim_mat(:)]
    dorsal_birthorder_similarity = [dorsal_birthorder_similarity; dorsal_pre(i).birth_order_sim(:)]
    dorsal_neurite_lengths = [dorsal_neurite_lengths; dorsal_pre(i).neurite_lengths]
end

self_index = dorsal_presynaptic_similarity == 1
dorsal_presynaptic_similarity(self_index) = []
dorsal_birthorder_similarity(self_index) = []
dorsal_neurite_lengths(self_index,:) = []

% Fit model
dorsal_presynaptic_data = table(dorsal_presynaptic_similarity, dorsal_birthorder_similarity, dorsal_neurite_lengths(:,1),dorsal_neurite_lengths(:,2),'VariableNames',{'Synapse_Similarity','Birthorder_Similarity','Youngest_Neuron','Mean_Age'})
lm_dorsal_pre = fitlm(dorsal_presynaptic_data,'Synapse_Similarity ~ Birthorder_Similarity + Mean_Age')
plot(lm_dorsal_pre)
ylim([0 1])
age_d = discretize(dorsal_presynaptic_data.Mean_Age,[0,11191.3082674945,17639.6673153414,21055.9963633213,27960.6969154466]);
map = [.25,.9,.9; .25,.25,.9; .9,.55,.25 ;.9,.25,.25];
%age_d = discretize(dorsal_presynaptic_data.Mean_Age,10)
%map = plasma(10)
figure; hold on
for i = 1:4
    scatter(dorsal_presynaptic_data.Birthorder_Similarity(age_d == i),dorsal_presynaptic_data.Synapse_Similarity(age_d==i),150,'MarkerFaceColor',map(i,:),'MarkerFaceAlpha',.75,'MarkerEdgeColor','k')
end

%% Birth Order Analysis of Dorsal Hemilineage Postsynapses

for i = 1:length(dorsal_hemilineages)
    if length(dorsal_hemilineages{i}) <= 2
    else
    temp_stats = birth_order_metricsv3(dorsal_hemilineages{i},2,0,savefigs,directory)
    dorsal_post(i).temporal_similarity = temp_stats.temporal_similarity;
    dorsal_post(i).birth_order_sim = temp_stats.birth_order_distance
    dorsal_post(i).sim_mat = temp_stats.sim_mat
    dorsal_post(i).dvs = temp_stats.dist_vs_sim
    dorsal_post(i).model = temp_stats.model
    dorsal_post(i).neurite_lengths = [temp_stats.min_distance_matrix(:),temp_stats.mean_distance_matrix(:)]

    clear temp_stats
    end 
end

% Compile synapse similarity, birthorder similarity, and average age
clear dorsal_postsynaptic_similarity and dorsal_birthorder_similarity and dorsal_neurite_lengths
dorsal_postsynaptic_similarity = dorsal_post(1).sim_mat(:)
dorsal_birthorder_similarity = dorsal_post(1).birth_order_sim(:)
dorsal_neurite_lengths = dorsal_post(1).neurite_lengths;
for i = 2:length(dorsal_post)
    dorsal_postsynaptic_similarity = [dorsal_postsynaptic_similarity; dorsal_post(i).sim_mat(:)]
    dorsal_birthorder_similarity = [dorsal_birthorder_similarity; dorsal_post(i).birth_order_sim(:)]
    dorsal_neurite_lengths = [dorsal_neurite_lengths; dorsal_post(i).neurite_lengths]
end

self_index = dorsal_postsynaptic_similarity == 1
dorsal_postsynaptic_similarity(self_index) = []
dorsal_birthorder_similarity(self_index) = []
dorsal_neurite_lengths(self_index,:) = []

% Fit model
dorsal_postsynaptic_data = table(dorsal_postsynaptic_similarity, dorsal_birthorder_similarity, dorsal_neurite_lengths(:,1),dorsal_neurite_lengths(:,2),'VariableNames',{'Synapse_Similarity','Birthorder_Similarity','Youngest_Neuron','Mean_Age'})
lm_dorsal_post = fitlm(dorsal_postsynaptic_data,'Synapse_Similarity ~ Birthorder_Similarity + Mean_Age')

map = [.25,.9,.9; .25,.25,.9; .9,.55,.25 ;.9,.25,.25];
figure; hold on
for i = 1:4
    scatter(dorsal_postsynaptic_data.Birthorder_Similarity(age_d == i),dorsal_postsynaptic_data.Synapse_Similarity(age_d==i),250,'MarkerFaceColor',map(i,:),'MarkerFaceAlpha',.75,'MarkerEdgeColor','k')
end


%% Birth Order Analysis of Ventral Hemilineage Presynapses

for i = 1:length(ventral_hemilineages)
    if length(ventral_hemilineages{i}) <= 4
    else
    [temp_stats] = birth_order_metricsv3(ventral_hemilineages{i},1,0,savefigs,directory)
    ventral_pre(i).temporal_similarity = temp_stats.temporal_similarity;
    ventral_pre(i).birth_order_sim = temp_stats.birth_order_distance;
    ventral_pre(i).sim_mat = temp_stats.sim_mat
    ventral_pre(i).dvs = temp_stats.dist_vs_sim
    ventral_pre(i).model = temp_stats.model
    ventral_pre(i).neurite_lengths = [temp_stats.min_distance_matrix(:),temp_stats.mean_distance_matrix(:)]

    clear temp_stats
    end
end

% Compile synapse similarity, birthorder similarity, and average age
clear ventral_presynaptic_similarity and dorsal_birthorder_similarity and dorsal_neurite_lengths
ventral_presynaptic_similarity = ventral_pre(1).sim_mat(:)
ventral_birthorder_similarity = ventral_pre(1).birth_order_sim(:)
ventral_neurite_lengths = ventral_pre(1).neurite_lengths;
for i = 2:length(ventral_pre)
    ventral_presynaptic_similarity = [ventral_presynaptic_similarity; ventral_pre(i).sim_mat(:)]
    ventral_birthorder_similarity = [ventral_birthorder_similarity; ventral_pre(i).birth_order_sim(:)]
    ventral_neurite_lengths = [ventral_neurite_lengths; ventral_pre(i).neurite_lengths]
end

self_index = ventral_presynaptic_similarity == 1
ventral_presynaptic_similarity(self_index) = []
ventral_birthorder_similarity(self_index) = []
ventral_neurite_lengths(self_index,:) = []

% Fit model
ventral_presynaptic_data = table(ventral_presynaptic_similarity, ventral_birthorder_similarity, ventral_neurite_lengths(:,1),ventral_neurite_lengths(:,2),'VariableNames',{'Synapse_Similarity','Birthorder_Similarity','Youngest_Neuron','Mean_Age'})
lm_ventral_pre = fitlm(ventral_presynaptic_data,'Synapse_Similarity ~ Birthorder_Similarity + Mean_Age')

age_d = discretize(ventral_presynaptic_data.Mean_Age,[0,11191.3082674945,17639.6673153414,21055.9963633213,27960.6969154466]);
map = [.25,.9,.9; .25,.25,.9; .9,.55,.25 ;.9,.25,.25];
figure; hold on
for i = 1:4
    scatter(ventral_presynaptic_data.Birthorder_Similarity(age_d == i),ventral_presynaptic_data.Synapse_Similarity(age_d==i),250,'MarkerFaceColor',map(i,:),'MarkerFaceAlpha',.75,'MarkerEdgeColor','k')
end

%% Birth Order Analysis of Ventral Hemilineage Postsynapses

for i = 1:length(ventral_hemilineages)
    if length(ventral_hemilineages{i}) <= 4
    else
    [temp_stats] = birth_order_metricsv3(ventral_hemilineages{i},2,0,savefigs,directory)
    ventral_post(i).temporal_similarity = temp_stats.temporal_similarity;
    ventral_post(i).birth_order_sim = temp_stats.birth_order_distance;
    ventral_post(i).sim_mat = temp_stats.sim_mat
    ventral_post(i).dvs = temp_stats.dist_vs_sim
    ventral_post(i).model = temp_stats.model
    ventral_post(i).neurite_lengths = [temp_stats.min_distance_matrix(:),temp_stats.mean_distance_matrix(:)]

    clear temp_stats
    end
end

% Compile synapse similarity, birthorder similarity, and average age
clear ventral_postsynaptic_similarity and ventral_birthorder_similarity and ventral_neurite_lengths
ventral_postsynaptic_similarity = ventral_post(1).sim_mat(:)
ventral_birthorder_similarity = ventral_post(1).birth_order_sim(:)
ventral_neurite_lengths = ventral_post(1).neurite_lengths;
for i = 2:length(ventral_post)
    ventral_postsynaptic_similarity = [ventral_postsynaptic_similarity; ventral_post(i).sim_mat(:)]
    ventral_birthorder_similarity = [ventral_birthorder_similarity; ventral_post(i).birth_order_sim(:)]
    ventral_neurite_lengths = [ventral_neurite_lengths; ventral_post(i).neurite_lengths]
end

self_index = ventral_postsynaptic_similarity == 1
ventral_postsynaptic_similarity(self_index) = []
ventral_birthorder_similarity(self_index) = []
ventral_neurite_lengths(self_index,:) = []



% Fit model
ventral_postsynaptic_data = table(ventral_postsynaptic_similarity, ventral_birthorder_similarity, ventral_neurite_lengths(:,1),ventral_neurite_lengths(:,2),'VariableNames',{'Synapse_Similarity','Birthorder_Similarity','Youngest_Neuron','Mean_Age'})
lm_ventral_post = fitlm(ventral_postsynaptic_data,'Synapse_Similarity ~ Birthorder_Similarity + Mean_Age')

map = [.25,.9,.9; .25,.25,.9; .9,.55,.25 ;.9,.25,.25];
figure; hold on
for i = 1:4
    scatter(ventral_postsynaptic_data.Birthorder_Similarity(age_d == i),ventral_postsynaptic_data.Synapse_Similarity(age_d==i),250,'MarkerFaceColor',map(i,:),'MarkerFaceAlpha',.75,'MarkerEdgeColor','k')
end

%%
allstat = [dorsal_pre,dorsal_post,ventral_pre,ventral_post];


for i = 1:length(allstat)
    if isempty(allstat(i).model);
    else
        
    coef(i) = allstat(i).model.Coefficients.pValue(2);
    rsq(i) = allstat(i).model.Rsquared.Adjusted;
    end
end
numel(find(coef<.05))/numel(coef)
mean(rsq)

figure; 
histogram(rsq,[-.2:.1:.6],'Normalization','cdf','FaceColor','k','FaceAlpha',.5)
xlabel('Per Lineage Birthorder vs Synapse Similarity R^2')
ylabel('Frequency')
set(gca,'FontSize',18)
%% Compile stats
for i = 1:length(dorsal_hemilineages)
    if isempty(dorsal_pre(i).model) == 0
        if dorsal_pre(i).model.Coefficients.pValue(2) < .05
            d_bom_pre{i} = dorsal_pre(i).birth_order_sim(:);
            d_som_pre{i} = dorsal_pre(i).sim_mat(:);
        else
        end
    else
    end
    if isempty(dorsal_post(i).model) == 0
        if dorsal_post(i).model.Coefficients.pValue(2) < .05
            d_bom_post{i} = dorsal_post(i).birth_order_sim(:);
            d_som_post{i} = dorsal_post(i).sim_mat(:);
        else
        end
    else
    end
end

for i = 1:length(ventral_hemilineages)
    if isempty(ventral_pre(i).model) == 0
        if ventral_pre(i).model.Coefficients.pValue(2) < .05
            v_bom_pre{i} = ventral_pre(i).birth_order_sim(:);
            v_som_pre{i} = ventral_pre(i).sim_mat(:);
        else
        end
    else
    end
    if isempty(ventral_post(i).model) == 0
        if ventral_post(i).model.Coefficients.pValue(2) < .05
            v_bom_post{i} = ventral_post(i).birth_order_sim(:);
            v_som_post{i} = ventral_post(i).sim_mat(:);
        else
        end
    else
    end
end

presom_v = cat(1,v_som_pre{:})
prebom_v = cat(1,v_bom_pre{:})
postsom_v = cat(1,v_som_post{:})
postbom_v = cat(1,v_bom_post{:})
clear v_bom_pre and v_som_pre and v_bom_post and v_som_post
 
presom_d = cat(1,d_som_pre{:})
prebom_d = cat(1,d_bom_pre{:})
postsom_d = cat(1,d_som_post{:})
postbom_d = cat(1,d_bom_post{:})

clear d_bom_pre and d_som_pre and d_bom_post and d_som_post
self = presom_v == 1;
presom_v(self) = []
prebom_v(self) = []
clear self
self = postsom_v == 1;
postsom_v(self) = []
postbom_v(self) = []
clear self

self = presom_d == 1;
presom_d(self) = []
prebom_d(self) = []
clear self
self = postsom_d == 1;
postsom_d(self) = []
postbom_d(self) = []
clear self

%% Total presynaptic birth order relationship
presynaptic_birth_order_metric = [prebom_d;prebom_v]
presynaptic_synapse_similarity = [presom_d;presom_v]
postsynaptic_birth_order_metric = [postbom_d;postbom_v]
postsynaptic_synapse_similarity = [postsom_d;postsom_v]
figure; hold on ; title('Brith Order Vs Presynaptic Similarity')
yyaxis left
scatter(presynaptic_birth_order_metric,presynaptic_synapse_similarity,200,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
f_total_pre = fitlm(presynaptic_birth_order_metric,presynaptic_synapse_similarity)
p = plot(f_total_pre)
ylabel('Presynaptic Similarity')
yyaxis right
scatter(postsynaptic_birth_order_metric,postsynaptic_synapse_similarity,200,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerEdgeAlpha',.8,'MarkerFaceAlpha',.5)
f_total_post = fitlm(postsynaptic_birth_order_metric,postsynaptic_synapse_similarity)
p2 = plot(f_total_post)
xlabel('Birth Order Similarity')
ylabel('Postsynaptic Similarity')
set(gca,'FontSize',18)
set(gcf,'Color','w')

% rmap = hsv(100)
% figure; hold on
% for i = 1:100
%     clear random_ind and randpresomb
%     random_ind = randperm(length(presomb))
%     randpresomb = presomb(random_ind)
%     scatter(prebomb,randpresomb,10,rmap(i,:))
% end


%% Total postsynaptic birth order relationship

figure; hold on ; title('Brith Order Vs Presynaptic Similarity')
scatter(presynaptic_birth_order_metric,presynaptic_synapse_similarity,100,'.','MarkerEdgeColor','r','MarkerFaceColor','k','MarkerEdgeAlpha',.8,'MarkerFaceAlpha',.5)
f_total_pre = fitlm(presynaptic_birth_order_metric,presynaptic_synapse_similarity)
p2 = plot(f_total_pre)
ylim([0,1])

xlabel('Birth Order Similarity')
ylabel('Presynaptic Similarity')
set(gca,'FontSize',18)
set(gcf,'Color','w')



figure; hold on ; title('Brith Order Vs Postsynaptic Similarity')
scatter(postsynaptic_birth_order_metric,postsynaptic_synapse_similarity,100,'.','MarkerEdgeColor','b','MarkerFaceColor','k','MarkerEdgeAlpha',.8,'MarkerFaceAlpha',.5)
f_total_post = fitlm(postsynaptic_birth_order_metric,postsynaptic_synapse_similarity)
p2 = plot(f_total_post)
ylim([0,1])

xlabel('Birth Order Similarity')
ylabel('Postsynaptic Similarity')
set(gca,'FontSize',18)
set(gcf,'Color','w')



all_synsimilarity = [dorsal_presynaptic_data.Synapse_Similarity; ventral_presynaptic_data.Synapse_Similarity]
all_birthsimilarity = [dorsal_presynaptic_data.Birthorder_Similarity;ventral_presynaptic_data.Birthorder_Similarity]
all_age = [dorsal_presynaptic_data.Mean_Age;ventral_presynaptic_data.Mean_Age]

all_age_d = discretize(all_age,[0,11191.3082674945,17639.6673153414,21055.9963633213,27960.6969154466])
map = [.25,.9,.9; .25,.25,.9; .9,.55,.25 ;.9,.25,.25];
figure; hold on
for i = 1:4
    scatter(all_birthsimilarity(all_age_d == i),all_synsimilarity(all_age_d==i),250,'MarkerFaceColor',map(i,:),'MarkerFaceAlpha',.75,'MarkerEdgeColor','k')
end


%% Put stats in an organized manner
% for i = 1:length(dorsal_hemilineages)
%     dnames{i} = [dorsal_hemilineages{i}(1).Annotations{:}]
% end
% col_names = {'Lineage','Pre','Post'}
% dorsal_model_statistics = table(transpose(dnames),transpose(dorsal_model_stats_pre),transpose(dorsal_model_stats_post),'VariableNames',col_names)
% 
% for i = 1:length(ventral_hemilineages)
%     vnames{i} = [ventral_hemilineages{i}(1).Annotations{:}]
% end
% 
% ventral_model_statistics = table(transpose(vnames),transpose(ventral_model_stats_pre),transpose(ventral_model_stats_post),'VariableNames',col_names)
% clear vnames and dnames and col_names
% % 
% % 
% % 

current_date = date
save(strcat(current_date,'Temporal_Targeting_Quantifications'))
