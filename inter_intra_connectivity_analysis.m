load mark2020_neurons_temporal_cohorts.mat
%%
% Index of temporal cohorts
temporal_cohort_options = unique(an_in(:,[2,4,5,6]))
[~,tc_index] = ismember(an_in(:,[2,4,5,6]),temporal_cohort_options,'rows')

% Index of hemilineages
hemilineage_options = unique(an_in(:,[2,4,5]))
[~,hl_index] = ismember(an_in(:,[2,4,5]),hemilineage_options,'rows')
%% Get Adjacency
adj = get_adjacency(nl,0,2)
%% Get a fractional degree matrix. That is, for each cohort, calculate the number of connections 
%  cohort i makes with cohort j / the number of possible connections between cohort i and cohort j

adj_tc = concatinate_adjacency(adj,nl,tc_index,1,1)
adj_hl = concatinate_adjacency(adj,nl,hl_index,1,1)
%% Get values of intra (diagonal) and inter (off diagonal) connection density

intra_tc = adj_tc(boolean(eye(size(adj_tc))));
inter_tc = adj_tc(boolean(-1*((eye(size(adj_tc)) - 1))));

intra_hl = adj_hl(boolean(eye(size(adj_hl))));
inter_hl = adj_hl(boolean(-1*((eye(size(adj_hl)) - 1))));


%% fraction of inter / intra connections

inter_tc_fraction = sum(double(inter_tc>0)) ./ numel(inter_tc)
inter_hl_fraction = sum(double(inter_hl>0)) ./ numel(inter_hl)

intra_tc_fraction = sum(double(intra_tc>0)) ./ numel(intra_tc)
intra_hl_fraction = sum(double(intra_hl>0)) ./ numel(intra_hl)


%% Plot

% Inter/Intra TC
figure; subplot(2,1,2); hold on
histogram(inter_tc(inter_tc>0),[0:.05:.5],'Normalization','probability','FaceColor','k')
histogram(intra_tc(intra_tc>0),[0:.05:.5],'Normalization','probability','FaceColor','r')

% Inter/Intra HL
subplot(2,1,1); hold on
histogram(inter_hl(inter_hl>0),[0:.1:1],'Normalization','probability','FaceColor','k')
histogram(intra_hl(intra_hl>0),[0:.1:1],'Normalization','probability','FaceColor','r')

% Inter TC/HL
figure;subplot(2,1,1); hold on
histogram(inter_tc(inter_tc>0),[0:.1:1],'Normalization','probability','FaceColor','c')
histogram(inter_hl(inter_hl>0),[0:.1:1],'Normalization','probability','FaceColor','r')
ylabel('Frequency')
xlabel('Interconnectivity Density')
set(gca,'FontSize',18)
title('')

% Intra TC/HL
subplot(2,1,2); hold on
histogram(intra_tc(intra_tc>0),[0:.1:1],'Normalization','probability','FaceColor','c')
histogram(intra_hl(intra_hl>0),[0:.1:1],'Normalization','probability','FaceColor','r')
ylabel('Frequency')
xlabel('Intraconnectivity Density')
set(gca,'FontSize',18)

%% Get the fractional combined adacency matrix:
% For each cohort i, how many synapses does it make with cohort j, divided
% by the total number of synapses recieved by cohort j. 
adj_tc_weighted = concatinate_adjacency(adj,nl,tc_index,0,1)
adj_hl_weighted = concatinate_adjacency(adj,nl,hl_index,0,1)
%% Get values of intra (diagonal) and inter (off diagonal) connection density

intra_tc_weighted = adj_tc_weighted(boolean(eye(size(adj_tc_weighted))));
inter_tc_weighted = adj_tc_weighted(boolean(-1*((eye(size(adj_tc_weighted)) - 1))));

intra_hl_weighted = adj_hl_weighted(boolean(eye(size(adj_hl_weighted))));
inter_hl_weighted = adj_hl_weighted(boolean(-1*((eye(size(adj_hl_weighted)) - 1))));

% Plot
figure;subplot(2,1,1); hold on
histogram(inter_tc_weighted(inter_tc_weighted>0),[0:.005:.1],'Normalization','probability','FaceColor','c')
histogram(inter_hl_weighted(inter_hl_weighted>0),[0:.005:.1],'Normalization','probability','FaceColor','r')
ylabel('Frequency')
xlabel('Interconnectivity Density')
set(gca,'FontSize',18)
title('')

subplot(2,1,2); hold on
histogram(intra_tc_weighted(intra_tc_weighted>0),[0:.005:.1],'Normalization','probability','FaceColor','c')
histogram(intra_hl_weighted(intra_hl_weighted>0),[0:.005:.1],'Normalization','probability','FaceColor','r')
ylabel('Frequency')
xlabel('Intraconnectivity Density')
set(gca,'FontSize',18)



