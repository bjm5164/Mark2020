load mark2020_neurons_temporal_cohorts.mat
% Remove NB1-2D from this analysis.
NB12D_index = boolean(strcmp(an_in.Lineage,'NB1-2') & an_in.DV_Index(:) == 1);
nl(NB12D_index) = [];
an_in(NB12D_index,:) = []; clear NB12D_index

% Get adjacency
adj = get_adjacency(nl,0,2)

% Remove self-self connections
adj(boolean(eye(size(adj)))) = 0

% Binarize
adj = double(adj>0);

%% Get temporal cohorts, hemilineages, and similarity

% Index of temporal cohorts
temporal_cohort_options = unique(an_in(:,[2,4,5,6]))
[~,tc_index] = ismember(an_in(:,[2,4,5,6]),temporal_cohort_options,'rows')

% Index of hemilineages
hemilineage_options = unique(an_in(:,[2,4,5]))
[~,hl_index] = ismember(an_in(:,[2,4,5]),hemilineage_options,'rows')

% Pre/Post similarity scores, thresholded. 
sim_thresh = [0:.1:.2]
for j = 1:length(sim_thresh);
    sim = synapse_similarity_io_overlap(nl,4000,4000);
    sim(boolean(eye(size(sim)))) = 0;
    sim(sim<sim_thresh(j)) = 0;
    similarity(:,:,j) = sim;

end

%% To look at inputs instead of outputs
adj = adj'
for i = 1:length(sim_thresh)
    similarity(:,:,i) = transpose(similarity(:,:,i));
end
%% Plot connection probability as a function of overlap. 

[N,edges,bin] = histcounts(similarity(:,:,1),[0:.1:1],'Normalization','probability')

for i = 1:10
    adj = adj;
    group = adj(bin == i)
    connectivity_fraction(i) = sum(group) ./ numel(group)
    group_count(i) = numel(group) / numel(similarity(:,:,1))
    number_connections(i) = sum(group) / sum(adj(:))
end

bins = [.05:.1:1]
figure; subplot(3,1,1)
bar(bins,group_count)
ylabel('Frequency')

subplot(3,1,2)
bar(bins,number_connections)
ylabel('Fraction of Connections')

subplot(3,1,3)
bar(bins,connectivity_fraction)
xlabel('Pre-Post Similarity')
ylabel('Fraction Connected')
    
    


%% Shuffle connections by weighted random sampling with weights being normalized pre-post overlap scores. 

for j = 1:length(sim_thresh)
 sim = similarity(:,:,j);   
    for i = 1:1000
        Wshuf_sim(:,:,i,j) = shufmat_sim(adj,sim);
        Wshuf_deg(:,:,i,j) = shufmat(adj);

    end
end

%% Graph statistics and global structure

% SVD for full real adj
explained_full_real = svd(adj);
var_full_real = (explained_full_real.^2)/sum(explained_full_real.^2);


% SVD for real temporal cohorts
combined_real = combined_adj(adj,tc_index);
%combined_real = double(combined_real>0);
explained_comb_real = svd(combined_real);
var_comb_real = (explained_comb_real.^2)/sum(explained_comb_real.^2);

% Node degrees for temporal cohorts

id_real = sum(combined_real,1);    % indegree = column sum of combadj
od_real = sum(combined_real,2)';   % outdegree = row sum of combadj



clear var_comb_shuf_deg and var_comb_shuf_sim and var_full_shuf_sim and var_full_shuf_deg

for j = 1:length(sim_thresh)
    for i = 1:1000
        % SVD on shuffuled full connectivity matrices with proximity based wiring probabilites 
        explained_fs = svd(Wshuf_sim(:,:,i,j));
        var_full_shuf_sim(:,i,j) = (explained_fs.^2)/sum(explained_fs.^2); 
        
        % SVD on binarized combined adjacency matrix with proximity based
        % wiring probabilities
        combined_shuf_sim = combined_adj(Wshuf_sim(:,:,i,j),tc_index);
        %combined_shuf_sim = double(combined_shuf_sim>0);
        explained_cs = svd(combined_shuf_sim);
        var_comb_shuf_sim(:,i,j) = (explained_cs.^2)/sum(explained_cs.^2);
        in_shuf = sum(combined_shuf_sim,1);    % indegree = column sum of combadj
        out_shuf = sum(combined_shuf_sim,2)';   % outdegree = row sum of combadj

        
        id_shuf_sim(:,i,j) = in_shuf; clear in_shuf
        od_shuf_sim(:,i,j) = out_shuf; clear out_shuf

        
        if j == 1
           
            explained_fd = svd(Wshuf_deg(:,:,i,j));
            var_full_shuf_deg(:,i) = (explained_fd.^2)/sum(explained_fd.^2);
            combined_shuf_deg = combined_adj(Wshuf_deg(:,:,i,j),tc_index);
            %combined_shuf_deg = double(combined_shuf_deg>0);
            explained_cd = svd(combined_shuf_deg);
            var_comb_shuf_deg(:,i) = (explained_cd.^2)/sum(explained_cd.^2);
            in_shuf = sum(combined_shuf_deg,1);
            out_shuf = sum(combined_shuf_deg,2);
            id_shuf_deg(:,i) = in_shuf; clear in_shuf
            od_shuf_deg(:,i) = out_shuf; clear out_shuf

        end
        
    end
end

%%

[n_id_real,edges] = histcounts(id_real,[0:1:20],'Normalization','probability')
[n_od_real,edges] = histcounts(od_real,[0:1:20],'Normalization','probability')

for i = 1:length(sim_thresh)
    for j = 1:1000
        n_in_shuf(:,j) = histcounts(id_shuf_sim(:,j,i),[0:1:20],'Normalization','probability');
        n_out_shuf(:,j) = histcounts(od_shuf_sim(:,j,i),[0:1:20],'Normalization','probability');
        if i == 1
            in_shuf_deg(:,j) = histcounts(id_shuf_deg,[0:1:20],'Normalization','probability');
            out_shuf_deg(:,j) = histcounts(od_shuf_deg,[0:1:20],'Normalization','probability');
        end
    end
    
    n_in_shuf_sim(:,i) = mean(n_in_shuf,2); 
    n_in_shuf_std(:,i) = std(n_in_shuf,[],2);  clear n_in_shuf
    
    n_out_shuf_sim(:,i) = mean(n_out_shuf,2); 
    n_out_shuf_std(:,i) = std(n_out_shuf,[],2); clear n_out_shuf
    
    if i == 1
        n_in_shuf_deg = mean(in_shuf_deg,2);
        std_in_shuf_deg = std(in_shuf_deg,[],2);
        
        n_out_shuf_deg = mean(out_shuf_deg,2);
        std_out_shuf_deg = std(out_shuf_deg,[],2);
    end
    
  
        
end
%%
map = cool(3)
bins = [.5:1:20]
figure; subplot(2,1,1); hold on
plot(bins,cumsum(n_od_real),'Color','r','LineWidth',3)
plot(bins,cumsum(n_out_shuf_deg),'Color','k')
errorbar(bins,cumsum(n_out_shuf_deg),std_in_shuf_deg,'Color','k','LineWidth',1)
for i = 1:3
    plot(bins,cumsum(n_out_shuf_sim(:,i)),'Color',map(i,:),'LineWidth',3)
    errorbar(bins,cumsum(n_out_shuf_sim(:,i)),n_out_shuf_std(:,i),'Color',map(i,:),'LineWidth',1)
end
ylabel('Frequency')
xlabel('Temporal Cohort Output Degree')
ylim([0 1])
set(gca,'FontSize',24)

subplot(2,1,2); hold on
plot(bins,cumsum(n_id_real),'Color','r','LineWidth',3)
plot(bins,cumsum(n_in_shuf_deg),'Color','k')
errorbar(bins,cumsum(n_in_shuf_deg),std_in_shuf_deg,'Color','k','LineWidth',1)
for i = 1:3
    plot(bins,cumsum(n_in_shuf_sim(:,i)),'Color',map(i,:),'LineWidth',3)
    errorbar(bins,cumsum(n_in_shuf_sim(:,i)),n_in_shuf_std(:,i),'Color',map(i,:),'LineWidth',1)
end
ylabel('Frequency')
xlabel('Temporal Cohort Input Degree')
ylim([0 1])

set(gca,'FontSize',24)


%%

figure;hold on
%plot([1:length(var_full_real)],var_full_real*100,'LineWidth',2,'Color','k');

x = [1:length(adj)]
for i = 1:length(sim_thresh)
    scatter(x,mean(var_full_shuf_sim(:,:,i)*100,2),100,'MarkerFaceColor',map(i,:),'MarkerEdgeColor',map(i,:),'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
    mean_shuf = mean(var_full_shuf_sim(:,:,i)*100,2);
    std_shuf = std(var_full_shuf_sim(:,:,i)*100,[],2);
    errorbar(x,mean_shuf,std_shuf,'Color',map(i,:),'LineWidth',3)
    
end

scatter(x,mean(var_full_shuf_deg*100,2),100,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5,'LineWidth',3)
mean_shuf_deg = mean(var_full_shuf_deg*100,2);
std_shuf_deg = std(var_full_shuf_deg*100,[],2)
for i = 1:length(x)
    errorbar(x(i),mean_shuf_deg(i),std_shuf_deg(i),'Color','k')
end

scatter([1:length(var_full_real)],var_full_real*100,100,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)

xlabel('Components')
ylabel('Percent Variance')
set(gca,'FontSize',18)
set(gcf,'Color','w')
%%
figure;hold on

x = [1:length(var_comb_real)]
for i = 1:length(sim_thresh)
    scatter(x,mean(var_comb_shuf_sim(:,:,i)*100,2),100,'MarkerFaceColor',map(i,:),'MarkerEdgeColor',map(i,:),'MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)
    mean_shuf_c = mean(var_comb_shuf_sim(:,:,i)*100,2);
    std_shuf_c = std(var_comb_shuf_sim(:,:,i)*100,[],2);
    for j = 1:length(x)
        errorbar(x(j),mean_shuf_c(j),std_shuf_c(j),'Color',map(i,:),'LineWidth',3)
    end
end

scatter(x,mean(var_comb_shuf_deg*100,2),100,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5,'LineWidth',3)
mean_shuf_deg_c = mean(var_comb_shuf_deg*100,2);
std_shuf_deg_c = std(var_comb_shuf_deg*100,[],2)
for i = 1:length(x)
    errorbar(x(i),mean_shuf_deg_c(i),std_shuf_deg_c(i),'Color','k')
end

scatter([1:length(var_comb_real)],var_comb_real*100,100,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerEdgeAlpha',.5,'MarkerFaceAlpha',.5)

xlabel('Components')
ylabel('Percent Variance')
set(gca,'FontSize',18)
set(gcf,'Color','w')
%% Cotargeting

real_co_target_prob = cohort_target_prob(adj,tc_index)

for i = 1:3
    shuf_cohort_target_prob(:,i) = cohort_target_prob(Wshuf_sim(:,:,:,i),tc_index)
end

shuf_cohort_target_prob_deg = cohort_target_prob(Wshuf_deg(:,:,:,1),tc_index)

%%
figure; hold on
for i = 1:3
    counts = histcounts(shuf_cohort_target_prob(:,i),[0:.02:.3],'Normalization','probability')
    area([.01:.02:.3],counts,'FaceColor',map(i,:),'FaceAlpha',.5); 
end

counts = histcounts(shuf_cohort_target_prob_deg,[0:.02:.3],'Normalization','probability')
area([.01:.02:.3],counts,'FaceColor','k','FaceAlpha',.5);

plot([real_co_target_prob, real_co_target_prob],[0, .35],'Color','r','LineWidth',3,'LineStyle','--')
xlabel('Co-target Probability')
ylabel('Frequency')
set(gca,'FontSize',18)
%%
newmap = {'m','b','c'}
figure;hold on
plot([1:length(var_comb_real)],var_comb_real*100,'LineWidth',2,'Color','k');

for i = 1:length(sim_thresh)
    
    shadedErrorBar([1:length(var_comb_shuf_sim(:,1,i))],mean(var_comb_shuf_sim(:,:,i)*100,2),std(var_comb_shuf_sim(:,:,i)*100,[],2),'lineprops',newmap{i})
end

shadedErrorBar([1:size(var_comb_shuf_deg,1)],mean(var_comb_shuf_deg*100,2),std(var_comb_shuf_deg*100,[],2),'lineprops','m')
% ylim([0 25])
% xlim([0 40])
xlabel('Components')
ylabel('Percent Variance')
set(gca,'FontSize',18)
set(gcf,'Color','w')
%% Look at degree distributions

real_out = sum(ad>0,2);
real_in = sum(adj>0,1);

for i = 1:3
    for j = 1:1000
        
        adj_shuf = Wshuf_sim(:,:,j,i);
        shuf_i(:,j,i) = sum(adj_shuf,2);
        shuf_o(:,j,i) = sum(adj_shuf,1);
    end
end

 
%%
for i = 1:3
    shuf_in(:,i) = mean(shuf_i(:,:,i),2);
    shuf_out(:,i) = mean(shuf_o(:,:,i),2);
end


%%
figure; hold on
for i = 1:5
    wf = Wshuf_sim(:,:,:,i);
    histogram(sum(wf,2))
end

%%
grid = [0 1 0 1 ; 1 0 1 0; 0 0 0 1; 1 1 1 0]
figure; subplot(1,2,1)
imagesc(grid); colormap(map)
axis off
grid_shuf = shufmat(grid)
subplot(1,2,2)
imagesc(grid_shuf); colormap(map)
axis off

sim = [.05 .09 .01 .25; .2 .1 .15 0; .05 .04 .01 .15; .2 .25 .15 .05]
map2 = plasma(50)
figure
imagesc(sim); colormap(cool)
grid_shuf_sim = shufmat_sim(grid,sim)
figure
imagesc(grid_shuf_sim); colormap(map)
%%
function Wshuf = shufmat(adj)
% Shuffle with input-degree based connection probabilities. 
    M = size(adj,1);
    N = size(adj,2);
    
    deg = sum(adj>0,2);
    
    cprobs = mean(adj > 0, 1);
    cprobs = cprobs / sum(cprobs);
    
    Wshuf = zeros(M,N);
    
    for i = 1:M

        if deg(i) > 0
            inds = datasample([1:N],deg(i),'Replace',false,'Weights',cprobs);
            Wshuf(i,inds) = 1;
        else
            Wshuf(i,:) = 0;
        end
    end

end


function Wshuf = shufmat_sim(adj,similarity)
% Shuffle with similarity score connection probabilities. 
    M = size(adj,1);
    N = size(adj,2);
    
    deg = sum(adj>0,2);
    
   
    Wshuf = zeros(M,N);
    
    for i = 1:M
        cprobs = similarity(i,:);
        cprobs = cprobs / sum(cprobs);
        if deg(i) > 0 & max(cprobs) > 0
            inds = datasample([1:N],deg(i),'Replace',false,'Weights',cprobs);
            Wshuf(i,inds) = 1;
        else
            Wshuf(i,:) = 0;
        end
    end

end


function comb_adj = combined_adj(adj,index)
% Combine adjacency matrices based on grouping index. 
    comb_adj = zeros(max(index));
    for i = 1:max(index);
        group_ind_i = find(index == i);
        for j = 1:max(index);
            group_ind_j = find(index == j);
            comb_adj(i,j) = sum(reshape(adj(group_ind_i,group_ind_j),numel(adj(group_ind_i,group_ind_j)),1));
            clear group_ind_j
        end
        clear group_ind_i
    end
end