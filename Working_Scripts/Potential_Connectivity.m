%% Analysis related to Peter's rule and temporal cohort connectivity specificity.
% First, we want to calculate both the input/output similarity of the
% neuron list as well as the potential and real adjacency matrices.

load mark2020_neurons_temporal_cohorts.mat
load 200507_potential_adjacency

% Potential adjacency w/ 5000um bandwidth for density calculation
potential_adj = potential_adjacency{5}
potential_adj(boolean(eye(size(potential_adj)))) = 0 % Set self-self connections to 0

% Similarity matrix
[sim_mat] = synapse_similarity_io_overlap(nl,2000)
sim_mat(boolean(eye(size(sim_mat)))) = 0 % Set self-self similarity to 0

% Fractional adjacency matrix
adj_frac = get_adjacency(nl,1)
adj_frac(adj_frac<.01) = 0;


% Real adjacency matrix 
adj_real = get_adjacency(nl,0)
adj_real(adj_real<2) = 0



%% Look at the number of connections

% Discretize similarity scores and count how many connections per exist between pairs with a given sim score. 
bins = [0:.2:1];
sim_mat_discrete = discretize(sim_mat,bins)
clear connections_per_bin and total_per_bin
for i = 1:length(bins)
    bin_idx = find(sim_mat_discrete == i);
    connections_per_bin(i) = nnz(adj_frac(bin_idx));
    total_per_bin(i) = numel(bin_idx)
end
connections_per_bin(isnan(connections_per_bin)) = 0

% Plot the connection density for each 
figure; subplot(2,1,1); hold on
area(bins,connections_per_bin./total_per_bin,'FaceColor','k','FaceAlpha',.5)
%xlim([0 .8])
set(gca,'FontSize',18)
xlabel('Pre-Post Similarity')
ylabel('Connection Density')

subplot(2,1,2)
H = histogram(sim_mat(:),bins,'FaceColor','k','FaceAlpha',.5)
E = H.BinEdges;
y = H.BinCounts;
xloc = E(1:end-1)+diff(E)/2;
text(xloc, y+1, string(y),'FontSize',18)
%xlim([0 .8])
ylim([0 3e4])
set(gca,'FontSize',18)
xlabel('Pre-Post Similarity')
ylabel('Frequency')

figure; hold on
scatter(sim_mat(:),adj_real(:))
lm = fitlm(sim_mat(:),adj_real(:))