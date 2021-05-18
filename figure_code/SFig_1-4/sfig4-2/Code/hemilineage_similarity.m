load mark2020_neurons_temporal_cohorts.mat
%%
ventral_index = unique([an_in.Lineage_Index(an_in.DV_Index == 0)])
ventral_legend = unique([an_in.Lineage(an_in.DV_Index == 0)])
ventral_hl = arrayfun(@(x) nl(an_in.Lineage_Index == ventral_index(x,1) & an_in.DV_Index == 0), 1:length(ventral_index),'UniformOutput',false)

dorsal_index = unique([an_in.Lineage_Index(an_in.DV_Index == 1)])
dorsal_legend = unique([an_in.Lineage(an_in.DV_Index == 1)])
dorsal_hl = arrayfun(@(x) nl(an_in.Lineage_Index == dorsal_index(x,1) & an_in.DV_Index == 1), 1:length(dorsal_index),'UniformOutput',false)
%%


save_path = '/Users/brandon/Documents/MATLAB/Repositories/Mark2020/sfig4-2/panels'

%%

for i = 1:max(an_in.Lineage_Index)
    clear dorsal and ventral
    if sum(dorsal_index == i) > 0
        dorsal = dorsal_hl{dorsal_index == i}
    else
        dorsal = []
    end
    if sum(ventral_index == i) > 0
        ventral = ventral_hl{ventral_index == i}
    else
        ventral = []
    end
    
    figure('pos',[1975 153 715 715],'rend','painters'); hold on
    
    if length(dorsal) > 0 | length(ventral) > 0
        if length(dorsal) == 0
            synapse_histograms(ventral,'b',2,.25)
            lin_title = unique(an_in.Lineage(an_in.Lineage_Index == i))
            saveas(gcf,string(strcat(save_path,'/',lin_title,'_postsynaptic','HL_Density.svg')),'svg')
        
        elseif length(ventral) == 0
            synapse_histograms(dorsal,'r',2,.25)
            lin_title = unique(an_in.Lineage(an_in.Lineage_Index == i))
            saveas(gcf,string(strcat(save_path,'/',lin_title,'_postsynaptic','HL_Density.svg')),'svg')
        else
            synapse_histograms(ventral,'b',2,.25)
            synapse_histograms(dorsal,'r',2,.25)
            lin_title = unique(an_in.Lineage(an_in.Lineage_Index == i))
            saveas(gcf,string(strcat(save_path,'/',lin_title,'_postsynaptic','HL_Density.svg')),'svg')

         
        end
    end
    
 
end
%%

% Parse lineages by lineage index and side index
unique_lineages = unique(an_in(:,[4]))
[~,lin_index] = ismember(an_in(:,[4]),unique_lineages(1:end,:),'rows')

lineages = arrayfun(@(x) nl(lin_index == x),1:max(lin_index),'UniformOutput',false)
lineage_legend = unique(an_in.Lineage)

lineages(7) = []
lineages(3) = []

lineage_legend(7) = []
lineage_legend(3) = []




for i = 1:length(lineages)
    lin_sim_pre{i} = synapse_similarity_v2(lineages{i},2000,3,[],1)
    lin_sim_post{i} = synapse_similarity_v2(lineages{i},2000,3,[],2);
end


%%
sim_to_compare = lin_sim_post;
dorsal_sim = []
ventral_sim = []
between_sim = []

for i = 1:length(lineages)
    lin_split = an_in(strcmp(an_in.Lineage,lineage_legend{i}),:);
    if length(unique(lin_split.DV_Index)) > 1
        
        dv_index = lin_split.DV_Index == 1
        lr_index = lin_split.Side_Index == 0

        dorsal_sim_l = reshape(sim_to_compare{i}(dv_index == 1 & lr_index == 1, dv_index == 1 & lr_index == 1),[],1);
        dorsal_sim_r = reshape(sim_to_compare{i}(dv_index == 1 & lr_index == 0, dv_index == 1 & lr_index == 0),[],1);
        dorsal_sim = [dorsal_sim;(dorsal_sim_l + dorsal_sim_r) * .5]
        dorsal_sim(dorsal_sim ==1) = []
        ds(i) = median(dorsal_sim)
        
        ventral_sim_l = reshape(sim_to_compare{i}(dv_index == 0 & lr_index == 1, dv_index == 0 & lr_index == 1),[],1);
        ventral_sim_r = reshape(sim_to_compare{i}(dv_index == 0 & lr_index == 0, dv_index == 0 & lr_index == 0),[],1);
        ventral_sim = [ventral_sim;(ventral_sim_l + ventral_sim_r) * .5]
        ventral_sim(ventral_sim == 1) = []
        vs(i) = median(ventral_sim)
        
        between_sim_l = reshape(sim_to_compare{i}(dv_index == 1 & lr_index == 1, dv_index == 0 & lr_index == 1),[],1);
        between_sim_r = reshape(sim_to_compare{i}(dv_index == 1 & lr_index == 0, dv_index == 0 & lr_index == 0),[],1);
        between_sim = [between_sim;(between_sim_l + between_sim_r) * .5];
        between_sim(between_sim == 1) = []
        bs(i) = median(between_sim)
    else
    end
end

%%
% figure; subplot(3,1,1)
% histogram(dorsal_sim,0:.1:1,'FaceColor','r','Normalization','pdf')
% subplot(3,1,2)
% histogram(ventral_sim,0:.1:1,'FaceColor', 'b','Normalization','pdf')
% subplot(3,1,3)
% histogram(between_sim,0:.1:1, 'FaceColor', 'k','Normalization','pdf')

figure('Position',[560 480 739 468]); hold on;
Violin(ds,1,'ViolinColor',[0 0 1],'Bandwidth',.1,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])
Violin(vs,2,'ViolinColor',[1 0,0],'Bandwidth',.1,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])
Violin(bs,3,'ViolinColor',[0 0 0],'Bandwidth',.1,'EdgeColor',[0,0,0],'BoxColor',[0,0,0])
xlim([.5,3.5])
xticks([1,2,3])
xticklabels({'Dorsal Hemilineages','Ventral Hemilineages','Between Hemilineages'})
ylabel('Presynapse Similarity')
set(gca,'FontSize',18)
%%
values = [ds,vs,bs]
groups = [ones(length(ds),1);ones(length(vs),1)+1;ones(length(bs),1)+2]
[p,tbl,stats] = anova1(values,groups)
multcompare(stats)

%%
nlist = lineages{1}(5:24)
nlist_l = nlist(1:2:end)
nlist_sim = synapse_similarity_v2(nlist_l,2000,3,[],1)

%%

for i = 1:length(lineages)
    nlist_sim_l_pre = synapse_similarity_v2(lineages{i}(1:2:end),2000,3,[],1)
    nlist_sim_r_pre = synapse_similarity_v2(lineages{i}(2:2:end),2000,3,[],1)
    nlist_sim_pre = (nlist_sim_l_pre + nlist_sim_r_pre)*.5
    
    nlist_sim_l_post = synapse_similarity_v2(lineages{i}(1:2:end),2000,3,[],2)
    nlist_sim_r_post = synapse_similarity_v2(lineages{i}(2:2:end),2000,3,[],2)
    nlist_sim_post = (nlist_sim_l_post + nlist_sim_r_post)*.5
    
    nlist_sim = (nlist_sim_pre + nlist_sim_post)*.5
    names = [lineages{i}.Names]
    
    Synapse_Distance_Clustering_v2(nlist_sim,names(1:2:end),2)
    saveas(gcf,string(strcat(save_path,'/',lineage_legend{i},'combined','_clustering.svg')),'svg')
end
    
