%% Parse lineages and hemilineages
load mark2020_neurons_temporal_cohorts.mat
[an_in,index] = sortrows(an_in,4);
nl = nl(index);

hemilineage_index = unique([an_in.Lineage_Index,an_in.DV_Index],'rows');

lineages = arrayfun(@(x) nl(an_in.Lineage_Index == x), 1:max(an_in.Lineage_Index),'UniformOutput',false);
hemilineages = arrayfun(@(x) nl(an_in.Lineage_Index == hemilineage_index(x,1) & an_in.DV_Index == hemilineage_index(x,2)), 1:length(hemilineage_index),'UniformOutput',false);

for i = 1:length(hemilineage_index)
    if hemilineage_index(i,2) == 0
      hemi_legend{i} = strcat(unique(an_in.Lineage(an_in.Lineage_Index == hemilineage_index(i,1))),'Ventral');
    elseif hemilineage_index(i,2) == 1
      hemi_legend{i} = strcat(unique(an_in.Lineage(an_in.Lineage_Index == hemilineage_index(i,1))),'Dorsal');
    else
      hemi_legend{i} = unique(an_in.Lineage(an_in.Lineage_Index == hemilineage_index(i,1)));
    end
end 

an_half = an_in(1:2:end,:);

savefigs = input('Save Figures? 1:Yes 0:no')
if savefigs == 1
    directory = uigetdir
else
    directory = 0
end

map = hsv(7)
map(end-1,:) = [1 0 1]
map(end-2,:) = [0.3 0 1.0000]
map(end-3,:) = [0 .6 1]
%% Calculate synapse similarity 

% Calculate presynapse similarity for neurons on left and right side.  
synaptic_similarity.presynaptic.left  = synapse_similarity_v2(nl(an_in.Side_Index == 0 & an_in.DV_Index < 2),2000,3,[],1);
synaptic_similarity.presynaptic.right = synapse_similarity_v2(nl(an_in.Side_Index == 1 & an_in.DV_Index < 2),2000,3,[],1);
synaptic_similarity.presynaptic.combined = (synaptic_similarity.presynaptic.left + synaptic_similarity.presynaptic.right)*.5 %Average similarities for left/right homologs


% Calculate postsynapse similarity for neurons on left and right sides.
synaptic_similarity.postsynaptic.left = synapse_similarity_v2(nl(an_in.Side_Index == 0 & an_in.DV_Index < 2),2000,3,[],2);
synaptic_similarity.postsynaptic.right = synapse_similarity_v2(nl(an_in.Side_Index == 1 & an_in.DV_Index < 2),2000,3,[],2);
synaptic_similarity.postsynaptic.combined = (synaptic_similarity.postsynaptic.left + synaptic_similarity.postsynaptic.right)*.5 % Average similarities for left/right homologs


% Get just dorsal similarities
synaptic_similarity.dorsal.presynaptic.left = synaptic_similarity.presynaptic.left(an_half.DV_Index == 1, an_half.DV_Index == 1);
synaptic_similarity.dorsal.presynaptic.right = synaptic_similarity.presynaptic.right(an_half.DV_Index == 1,an_half.DV_Index == 1);
synaptic_similarity.dorsal.presynaptic.combined = (synaptic_similarity.dorsal.presynaptic.left + synaptic_similarity.dorsal.presynaptic.right)*.5;

synaptic_similarity.dorsal.postsynaptic.left = synaptic_similarity.postsynaptic.left(an_half.DV_Index == 1,an_half.DV_Index == 1);
synaptic_similarity.dorsal.postsynaptic.right = synaptic_similarity.postsynaptic.right(an_half.DV_Index == 1,an_half.DV_Index == 1);
synaptic_similarity.dorsal.postsynaptic.combined = (synaptic_similarity.dorsal.postsynaptic.left + synaptic_similarity.dorsal.postsynaptic.right)*.5;

synaptic_similarity.dorsal.combined = (synaptic_similarity.dorsal.presynaptic.combined+synaptic_similarity.dorsal.postsynaptic.combined)*.5;

% Get just ventral similarities
synaptic_similarity.ventral.presynaptic.left = synaptic_similarity.presynaptic.left(an_half.DV_Index == 0,an_half.DV_Index == 0);
synaptic_similarity.ventral.presynaptic.right = synaptic_similarity.presynaptic.right(an_half.DV_Index == 0,an_half.DV_Index == 0);
synaptic_similarity.ventral.presynaptic.combined = (synaptic_similarity.ventral.presynaptic.left + synaptic_similarity.ventral.presynaptic.right)*.5

synaptic_similarity.ventral.postsynaptic.left = synaptic_similarity.postsynaptic.left(an_half.DV_Index == 0,an_half.DV_Index == 0);
synaptic_similarity.ventral.postsynaptic.right = synaptic_similarity.postsynaptic.right(an_half.DV_Index == 0,an_half.DV_Index == 0);
synaptic_similarity.ventral.postsynaptic.combined = (synaptic_similarity.ventral.postsynaptic.left + synaptic_similarity.ventral.postsynaptic.right)*.5

synaptic_similarity.ventral.combined = (synaptic_similarity.ventral.presynaptic.combined+synaptic_similarity.ventral.postsynaptic.combined)*.5;
%%

% Average similarity for pre and post synapses
synaptic_similarity.combined_similarity = (synaptic_similarity.presynaptic.combined+synaptic_similarity.postsynaptic.combined)*.5

for i = 1:length(synaptic_similarity.presynaptic.left)
    precorco = corrcoef(synaptic_similarity.presynaptic.left(i,:),synaptic_similarity.presynaptic.right(i,:))
    pre_cor(i) = precorco(1,2)
    
    postcorco = corrcoef(synaptic_similarity.postsynaptic.left(i,:),synaptic_similarity.postsynaptic.right(i,:))
    post_cor(i) = postcorco(1,2)
    
    lr_diff_pre(i,:) = synaptic_similarity.presynaptic.left(i,:)-synaptic_similarity.presynaptic.right(i,:)
    lr_diff_post(i,:) = synaptic_similarity.postsynaptic.left(i,:)-synaptic_similarity.postsynaptic.right(i,:)
end

figure; hold on
histogram(pre_cor,.7:.01:1,'FaceColor','c','Normalization','probability')
histogram(post_cor,.7:.01:1,'FaceColor','r','Normalization','probability')
xlabel('Left / Right Homolog correlation coefficient')
ylabel('Frequency')
legend({'Presynaptic Similarity','Postsynaptic Similarity'})
set(gca,'FontSize',18)
[h,p] = kstest2(pre_cor,post_cor)


figure; hold on 
histogram(lr_diff_pre(:),-.45:.02:.45,'FaceColor','c','Normalization','probability')
histogram(lr_diff_post(:),-.45:.02:.45,'FaceColor','r','Normalization','probability')
xlabel('Left / Right Homolog Difference')
ylabel('Frequency')
legend({'Presynaptic Similarity','Postsynaptic Similarity'})
set(gca,'FontSize',18)
[h,p] = kstest2(lr_diff_pre(:),lr_diff_post(:))

fraction_larger_pre = sum(abs(lr_diff_pre(:))>.03)/numel(lr_diff_pre)
fraction_larger_post = sum(abs(lr_diff_post(:))>.01)/numel(lr_diff_post)

%% Calculate covarience matrix for left/right 
cov_mat = corrcov(synaptic_similarity.presynaptic.left,synaptic_similarity.postsynaptic.right)

%%
% Remove self-comparisons
synaptic_similarity.presynaptic.left(synaptic_similarity.presynaptic.left == 1) = []
synaptic_similarity.presynaptic.right(synaptic_similarity.presynaptic.right == 1) = []
synaptic_similarity.postsynaptic.left(synaptic_similarity.postsynaptic.left == 1) = []
synaptic_similarity.postsynaptic.right(synaptic_similarity.postsynaptic.right == 1) = []

synaptic_similarity.dorsal.presynaptic.left(synaptic_similarity.dorsal.presynaptic.left == 1) = []
synaptic_similarity.dorsal.presynaptic.right(synaptic_similarity.dorsal.presynaptic.right == 1) = []
synaptic_similarity.dorsal.postsynaptic.left(synaptic_similarity.dorsal.postsynaptic.left == 1) = []
synaptic_similarity.dorsal.postsynaptic.right(synaptic_similarity.dorsal.postsynaptic.right == 1) = []

synaptic_similarity.ventral.presynaptic.left(synaptic_similarity.ventral.presynaptic.left == 1) = []
synaptic_similarity.ventral.presynaptic.right(synaptic_similarity.ventral.presynaptic.right == 1) = []
synaptic_similarity.ventral.postsynaptic.left(synaptic_similarity.ventral.postsynaptic.left == 1) = []
synaptic_similarity.ventral.postsynaptic.right(synaptic_similarity.ventral.postsynaptic.right == 1) = []

% Calculate left/right correlation coefficients. 
pre_cor = corrcoef(synaptic_similarity.presynaptic.left,synaptic_similarity.presynaptic.right)
post_cor = corrcoef(synaptic_similarity.postsynaptic.left,synaptic_similarity.postsynaptic.right)




% Plot left vs right pre/post synapse similarity.
figure; subplot(2,1,1)
scatter(synaptic_similarity.presynaptic.left(:),synaptic_similarity.presynaptic.right(:),100,'o','MarkerFaceAlpha',.15,'MarkerEdgeAlpha',.25,'MarkerFaceColor','r','MarkerEdgeColor','k')
xlabel('Synapse Similarity Left')
ylabel('Synapse Similarity Right')
set(gca,'FontSize',14)
axis equal
xlim([0,1])
ylim([0,1])
legend({'Presynapses'})
subplot(2,1,2)
scatter(synaptic_similarity.postsynaptic.left(:),synaptic_similarity.postsynaptic.right(:),100,'o','MarkerFaceAlpha',.15,'MarkerEdgeAlpha',.25,'MarkerFaceColor','c','MarkerEdgeColor','k')
xlabel('Synapse Similarity Left')
ylabel('Synapse Similarity Right')
set(gca,'FontSize',14)
axis equal
xlim([0,1])
ylim([0,1])
legend({'Postsynapses'})
%%
clear pre_sim_l and pre_sim_r and pre_sim_l_d and pre_sim_r_d and pre_sim_l_v and pre_sim_r_v and post_sim_l and post_sim_r and post_sim_l_d ...
    and post_sim_r_d and post_sim_l_v and post_sim_r_v 
%% Look at pre- and post- synaptic similarity of dorsal and ventral hemilineages 
dorsal_neurons = nl(an_in.DV_Index == 1)
dorsal_an = an_half(an_half.DV_Index == 1,:)
dorsal_index = unique(dorsal_an(:,[2,4]),'rows');

% Get presynapse similarity for dorsal hemilineages
similarity_matrix = synaptic_similarity.dorsal.presynaptic.combined
for i = 1:length(dorsal_index.DV_Index)
    related_index_h = ismember(dorsal_an(:,[2,4]),dorsal_index(i,:),'rows'); % Find neurons in hemilineage i
    unrelated_index_h = boolean((ismember(dorsal_an(:,4),dorsal_index(i,2))*-1)+1)
    sr_h = similarity_matrix(related_index_h,related_index_h); % get pairwise similarities between them
    sur_h = similarity_matrix(related_index_h,unrelated_index_h);
    sim_related = sr_h(:); clear sr_h
    sim_unrelated = sur_h(:); clear sur_h
    sim_related(sim_related == 1) = []; % remove self-self comparisons
    dorsal_presynaptic_similarity{i} = sim_related(:); 
    dorsal_presynaptic_similarity_unrelated{i} = sim_unrelated(:); clear sim_related and related_index_h and sim_unrelated and unrelated_index_h
end
clear similarity_matrix 

% Get postsynapse similarity for dorsal hemilineages
similarity_matrix = synaptic_similarity.dorsal.postsynaptic.combined
for i = 1:length(dorsal_index.DV_Index)
    related_index_h = ismember(dorsal_an(:,[2,4]),dorsal_index(i,:),'rows'); % Find neurons in hemilineage i
    unrelated_index_h = boolean((ismember(dorsal_an(:,4),dorsal_index(i,2))*-1)+1)
    sr_h = similarity_matrix(related_index_h,related_index_h); % get pairwise similarities between them
    sur_h = similarity_matrix(related_index_h,unrelated_index_h)
    sim_related = sr_h(:); clear sr_h
    sim_unrelated = sur_h(:); clear sur_h
    sim_related(sim_related == 1) = []; % remove self-self comparisons
    dorsal_postsynaptic_similarity{i} = sim_related(:); 
    dorsal_postsynaptic_similarity_unrelated{i} = sim_unrelated(:); clear sim_related and related_index_h and sim_unrelated and unrelated_index_h
end
clear similarity_matrix

ventral_neurons = nl(an_in.DV_Index == 0);
ventral_an = an_half(an_half.DV_Index == 0,:);
ventral_index = unique(ventral_an(:,[2,4]),'rows'); 

% Get presynapse similarity for ventral hemilineages
similarity_matrix = synaptic_similarity.ventral.presynaptic.combined
for i = 1:length(ventral_index.DV_Index)
    related_index_h = ismember(ventral_an(:,[2,4]),ventral_index(i,:),'rows'); % Find neurons in hemilineage i
    unrelated_index_h = boolean((ismember(ventral_an(:,4),ventral_index(i,2))*-1)+1);
    sr_h = similarity_matrix(related_index_h,related_index_h); % get pairwise similarities between them
    sur_h = similarity_matrix(related_index_h,unrelated_index_h);
    sim_related = sr_h(:); clear sr_h
    sim_unrelated = sur_h(:); clear sur_h
    sim_related(sim_related == 1) = []; % remove self-self comparisons
    ventral_presynaptic_similarity{i} = sim_related(:); 
    ventral_presynaptic_similarity_unrelated{i} = sim_unrelated(:);clear sim_related and related_index_h and sim_unrelated and unrelated_index_h
end
clear similarity_matrix

% Get postsynapse similarity for ventral hemilineages 
similarity_matrix = synaptic_similarity.ventral.postsynaptic.combined
for i = 1:length(ventral_index.DV_Index)
    related_index_h = ismember(ventral_an(:,[2,4]),ventral_index(i,:),'rows'); % Find neurons in hemilineage i
    unrelated_index_h = boolean((ismember(ventral_an(:,4),ventral_index(i,2))*-1)+1);
    sr_h = similarity_matrix(related_index_h,related_index_h); % get pairwise similarities between them
    sur_h = similarity_matrix(related_index_h,unrelated_index_h);
    sim_related = sr_h(:); clear sr_h
    sim_unrelated = sur_h(:); clear sur_h
    sim_related(sim_related == 1) = []; % remove self-self comparisons
    ventral_postsynaptic_similarity{i} = sim_related(:); 
    ventral_postsynaptic_similarity_unrelated{i} = sim_unrelated(:);clear sim_related and related_index_h and sim_unrelated and unrelated_index_h
end
clear similarity_matrix

dorsal_presynaptic = cat(1,dorsal_presynaptic_similarity{:});
dorsal_presynaptic_unrelated = cat(1,dorsal_presynaptic_similarity_unrelated{:});

dorsal_postsynaptic = cat(1,dorsal_postsynaptic_similarity{:});
dorsal_postsynaptic_unrelated = cat(1,dorsal_postsynaptic_similarity_unrelated{:});

ventral_presynaptic = cat(1,ventral_presynaptic_similarity{:})
ventral_presynaptic_unrelated = cat(1,ventral_presynaptic_similarity_unrelated{:});

ventral_postsynaptic = cat(1,ventral_postsynaptic_similarity{:})
ventral_postsynaptic_unrelated = cat(1,ventral_postsynaptic_similarity_unrelated{:})

%%
figure; hold on
bar(1,nanmean(dorsal_presynaptic_unrelated),'FaceColor','k')
errorbar(1,nanmean(dorsal_presynaptic_unrelated),nanstd(dorsal_presynaptic_unrelated)/sqrt(numel(dorsal_presynaptic_unrelated)),'k')
bar(2,nanmean(dorsal_presynaptic),'FaceColor','b')
errorbar(2,nanmean(dorsal_presynaptic),nanstd(dorsal_presynaptic)/sqrt(numel(dorsal_presynaptic)),'k')
[p_d_pre,h] = ranksum(dorsal_presynaptic,dorsal_presynaptic_unrelated)

bar(7,nanmean(dorsal_postsynaptic_unrelated),'FaceColor','k')
errorbar(7,nanmean(dorsal_postsynaptic_unrelated),std(dorsal_postsynaptic_unrelated)/sqrt(numel(dorsal_postsynaptic_unrelated)),'k')
bar(8,nanmean(dorsal_postsynaptic),'FaceColor',[0.9255    0.3765    0.0471])
errorbar(8,nanmean(dorsal_postsynaptic),nanstd(dorsal_postsynaptic)/sqrt(numel(dorsal_postsynaptic)),'k')
[p_d_post,h] = ranksum(dorsal_postsynaptic,dorsal_postsynaptic_unrelated)

bar(4,nanmean(ventral_presynaptic_unrelated),'FaceColor','k')
errorbar(4,nanmean(ventral_presynaptic_unrelated),std(ventral_presynaptic_unrelated)/sqrt(numel(ventral_presynaptic_unrelated)),'k')
bar(5,nanmean(ventral_presynaptic),'FaceColor','b')
errorbar(5,nanmean(ventral_presynaptic),nanstd(ventral_presynaptic)/sqrt(numel(ventral_presynaptic)),'k')
[p_v_pre,h] = ranksum(ventral_presynaptic,ventral_presynaptic_unrelated)

bar(10,nanmean(ventral_postsynaptic_unrelated),'FaceColor','k')
errorbar(10,nanmean(ventral_postsynaptic_unrelated),std(ventral_postsynaptic_unrelated)/sqrt(numel(ventral_postsynaptic_unrelated)),'k')
bar(11,nanmean(ventral_postsynaptic),'FaceColor',[0.9255    0.3765    0.0471])
errorbar(11,nanmean(ventral_postsynaptic),std(ventral_postsynaptic)/sqrt(numel(ventral_postsynaptic)),'k')
[p_v_post,h] = ranksum(ventral_postsynaptic,ventral_postsynaptic_unrelated)

xticks([1.5,4.5,7.5,10.5])
xticklabels({'Dorsal Presynaptic','Ventral Presynaptic','Dorsal Postsynaptic','Ventral Postsynaptic'})
ylabel('Synapse Similarity')
set(gca,'FontSize',18)

display('dorsal pre')
nanmean(dorsal_presynaptic)/nanmean(dorsal_presynaptic_unrelated)
display('dorsal post')
nanmean(dorsal_postsynaptic)/nanmean(dorsal_postsynaptic_unrelated)
display('ventral pre')
nanmean(ventral_presynaptic)/nanmean(ventral_presynaptic_unrelated)
display('ventral_post')
nanmean(ventral_postsynaptic)/nanmean(ventral_postsynaptic_unrelated)


%% Cluster Dorsal Hemilineages
lindex_d = an_in.Lineage_Index(an_in.DV_Index==1 & an_in.Side_Index == 1) % Get lineage index for dorsal hemilineages on one side.
d_lin = unique(lindex_d) 

[h1 t1 permd Z] = Synapse_Distance_Clustering_v2(synaptic_similarity.dorsal.combined,an_in.Names(an_in.DV_Index==1 & an_in.Side_Index == 1 ),7,0) % Cluster dorsal neurons by synapse similarity.
reordered_hemi = lindex_d(permd);

% Recolor labels to match lineages
subplot(3,4,[1 5 9]);
ax = get(gca);
lab = ax.YAxis.TickLabels;
for i = 1:length(lab)
    
    end_name_ind = strfind(lab{i},',');
    lab{i}(end_name_ind:end) = [];
end

loc = ax.YAxis.TickValues;
for k = 1:numel(lab) % for every type of lable
    ind = strcmp(lab{k},lab); % find all instances in lab
    x = repelem(ax.XAxis.Limits(1)-0.01,sum(ind)); % make an x position vector
     % place this lable at the same locations with a distinct color:
    text(x,loc(ind),lab(ind),'Color',map(reordered_hemi(k)-1,:));
    ax.YAxis.TickLabels = []; % remove the original labels
% replace the original labels with white space, to keep the axes position:
ax.YAxis.TickLabels = repelem('  ',max(cellfun(@numel,lab)));
set(gca,'FontSize',18)
end

if savefigs == 1
    saveas(gcf,strcat(directory,'/','Dorsal_Hemilineage_Similarity_Combined'),'svg')
    close all
else
end

%% Cluster Ventral Hemilineages
lindex_v = an_in.Lineage_Index(an_in.DV_Index==0 & an_in.Side_Index == 1) % Get lineage index for ventral hemilienages on one side.
v_lin = unique(lindex_v)

[h1 t1 permv Z] = Synapse_Distance_Clustering_v2(synaptic_similarity.ventral.combined,an_in.Names(an_in.DV_Index==0 & an_in.Side_Index == 1 ),7,0) % Cluster ventral neurons by synapse similarity.
reordered_hemi = lindex_v(permv);

% Recolor labels to match lineages.
subplot(3,4,[1 5 9]);
ax = get(gca);
lab = ax.YAxis.TickLabels;
for i = 1:length(lab)
    
    end_name_ind = strfind(lab{i},',');
    lab{i}(end_name_ind:end) = [];
end

loc = ax.YAxis.TickValues;
for k = 1:numel(lab) % for every type of lable
    ind = strcmp(lab{k},lab); % find all instances in lab
    x = repelem(ax.XAxis.Limits(1)-0.01,sum(ind)); % make an x position vector
     % place this lable at the same locations with a distinct color:
    text(x,loc(ind),lab(ind),'Color',map(reordered_hemi(k)-1,:));
    ax.YAxis.TickLabels = []; % remove the original labels
% replace the original labels with white space, to keep the axes position:
ax.YAxis.TickLabels = repelem('  ',max(cellfun(@numel,lab)));
set(gca,'FontSize',18)
end


if savefigs == 1
    saveas(gcf,strcat(directory,'/','Ventral_Hemilineage_Similarity_Combined'),'svg')
    close all
else
end


%% Similarity Matrices 
% Generates one half at a time with different colormaps.  Switch lines
% to generate upper and lower halves for pre/post

%Dorsal pre-similarity
lindex_d = an_in.Lineage_Index(an_in.DV_Index==1 & an_in.Side_Index == 1) % Get lineage index for dorsal hemilineages on one side.

figure('pos',[1400,100,1200,500]); 
subplot(1,2,1)
M = ones(size(synaptic_similarity.dorsal.presynaptic.combined))
M_up = boolean(tril(M))
M_low = boolean(triu(M))
sim_d = nan(size(synaptic_similarity.dorsal.presynaptic.combined))
%sim_d(M_up) = synaptic_similarity.dorsal.presynaptic.combined(M_up) %PRE
sim_d(M_low) = synaptic_similarity.dorsal.postsynaptic.combined(M_low); %POST
clear M and M_up and M_low

imagesc(sim_d,[0,.5]); axis xy
yticks(1:length(synaptic_similarity.dorsal.presynaptic.combined))
yticklabels(an_in.Names(an_in.DV_Index==1 & an_in.Side_Index == 1 ))
xticks([])
ax = get(gca);
lab = ax.YAxis.TickLabels;
for i = 1:length(lab)
    
    end_name_ind = strfind(lab{i},',');
    lab{i}(end_name_ind:end) = [];
end

loc = ax.YAxis.TickValues;
for k = 1:numel(lab) % for every type of lable
    ind = strcmp(lab{k},lab); % find all instances in lab
    x = repelem(ax.XAxis.Limits(1)-2.5,sum(ind)); % make an x position vector
     % place this lable at the same locations with a distinct color:
    text(x,loc(ind),lab(ind),'Color',map(lindex_d(k)-1,:));
    ax.YAxis.TickLabels = []; % remove the original labels
% replace the original labels with white space, to keep the axes position:
ax.YAxis.TickLabels = repelem('  ',max(cellfun(@numel,lab)));
set(gca,'FontSize',18)
end
%colormap(cbrewer('seq','Blues',64)) % PRE
colormap(cbrewer('seq','Oranges',64)) % POST

title('Dorsal Synaptic Similarity')


% Ventral pre-similarity
lindex_v = an_in.Lineage_Index(an_in.DV_Index==0 & an_in.Side_Index == 1) % Get lineage index for ventral hemilienages on one side.
subplot(1,2,2)
M = ones(size(synaptic_similarity.ventral.presynaptic.combined))
M_up = boolean(tril(M))
M_low = boolean(triu(M))
sim_v = nan(size(synaptic_similarity.ventral.presynaptic.combined)) 
%sim_v(M_up) = synaptic_similarity.ventral.presynaptic.combined(M_up) %PRE
sim_v(M_low) = synaptic_similarity.ventral.postsynaptic.combined(M_low); %POST
clear M and M_up and M_low 
imagesc(sim_v,[0,.5]); axis xy
axis xy
yticks(1:length(synaptic_similarity.ventral.presynaptic.combined))
yticklabels(an_in.Names(an_in.DV_Index==0 & an_in.Side_Index == 1 ))
xticks([])
ax = get(gca);
lab = ax.YAxis.TickLabels;
for i = 1:length(lab)
    
    end_name_ind = strfind(lab{i},',');
    lab{i}(end_name_ind:end) = [];
end

loc = ax.YAxis.TickValues;
for k = 1:numel(lab) % for every type of lable
    ind = strcmp(lab{k},lab); % find all instances in lab
    x = repelem(ax.XAxis.Limits(1)-2.5,sum(ind)); % make an x position vector
     % place this lable at the same locations with a distinct color:
    text(x,loc(ind),lab(ind),'Color',map(lindex_v(k)-1,:));
    ax.YAxis.TickLabels = []; % remove the original labels
% replace the original labels with white space, to keep the axes position:
ax.YAxis.TickLabels = repelem('  ',max(cellfun(@numel,lab)));
set(gca,'FontSize',18)
end
c = colorbar
c.Label.String = 'Synapse Similarity'
%colormap(cbrewer('seq','Blues',64)) % PRE
colormap(cbrewer('seq','Oranges',64)) % POST
title('Ventral Synaptic Similarity')



