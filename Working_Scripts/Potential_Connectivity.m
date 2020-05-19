%% Analysis related to Peter's rule and temporal cohort connectivity specificity.
% First, we want to calculate both the input/output similarity of the
% neuron list as well as the potential and real adjacency matrices.

load mark2020_neurons_temporal_cohorts.mat
load 200507_potential_adjacency

% Potential adjacency w/ 5000um bandwidth for density calculation
potential_adj = potential_adjacency{5}
potential_adj(boolean(eye(size(potential_adj)))) = 0 % Set self-self connections to 0

% Similarity matrix
[sim_mat] = synapse_similarity_io_overlap(nl,2000,5000)
sim_pre_to_post = sim_mat; % sim_mat(i,j) is how similar presynapses of neuron i are to postsynapses of neuron j
sim_post_to_pre = sim_mat'; % sim_mat(i,j) is how similar postsynapses of neuron i are to presynapses of neuron j
sim_mat(boolean(eye(size(sim_mat)))) = 0 % Set self-self similarity to 0

% Fractional adjacency matrix
adj_frac = get_adjacency(nl,1)
adj_frac(adj_frac<.01) = 0;


% Real adjacency matrix 
adj_real = get_adjacency(nl,0)
adj_real(adj_real<2) = 0

%% Look at number of synapses versus similarity.
% The number of zeros is hard to read on a scatter, so let's remove them
% and plot them separately as a histogram.

connected_idx = find(adj_frac(:)>0);

figure; subplot(2,1,1); hold on
scatter(sim_mat(connected_idx),adj_frac(connected_idx),20,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',.5)
xlabel('Pre/Post Similarity')
ylabel('Connection Strength (fraction of total iput)')
%lm = fitlm(sim_mat(connected_idx),adj_real(connected_idx))
%plot(lm)
set(gca,'FontSize',14)
xlim([0,1])
subplot(2,1,2);
histogram(sim_mat(setdiff([1:numel(sim_mat)],connected_idx)),0:.1:1,'FaceColor','k','FaceAlpha',.5,'Normalization','probability')
xlabel('Pre/Post Similarity')
ylabel('Frequency')
set(gca,'FontSize',14)
xlim([0,1])
set(gca,'YScale','log')



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
area(bins+.1,connections_per_bin./total_per_bin,'FaceColor','k','FaceAlpha',.5)
xlim([0 1])
xticks([0:.2:1])
set(gca,'FontSize',18)
xlabel('Pre-Post Similarity')
ylabel('Connection Density')

subplot(2,1,2)
H = histogram(sim_mat(:),bins,'FaceColor','k','FaceAlpha',.5)
E = H.BinEdges;
y = H.BinCounts;
xloc = E(1:end-1)+diff(E)/2;
text(xloc, y+100, string(y),'FontSize',18)
xlim([0 1])
%ylim([0 3e4])
set(gca,'FontSize',18)
xlabel('Pre-Post Similarity')
ylabel('Frequency')
set(gca,'YScale','log')

%% Now look at connection probabilities for temporal cohorts.

neuron_deg = adj_frac
neuron_deg(neuron_deg>0) = 1
sim_to_use = sim_pre_to_post;
sim_thresh = [0:.1:.3]

[co_target_prob,co_target_fraction,neurons_used,similarity_distributions] = temporal_cohort_co_targeting(nl,an_in,sim_to_use,neuron_deg)
[co_target_prob_rand,co_target_fraction_rand,neurons_used_rand,similarity_distributions_random] = temporal_cohort_co_targeting(nl,an_in,sim_to_use,neuron_deg,1,sim_thresh)


%%

map = cbrewer('seq','Blues',length(co_target_prob_rand))

figure; hold on
for i = 1:length(co_target_prob_rand)
    random_counts = histcounts(sum(co_target_prob_rand{i},2)/100,0:.04:1,'Normalization','probability')
    area([0:.04:.98],random_counts,'FaceColor',map(i,:),'FaceAlpha',.75,'LineWidth',1.5,'EdgeColor',map(i,:)*.5)
    %histogram(sum(same_cohort_random{i},2)/100,0:.02:1,'FaceColor',map(i,:),'Normalization','probability')
end

real_counts = histcounts(sum(co_target_prob,2)/100,0:.04:1,'Normalization','probability')
plot([0:.04:.98],real_counts,'Color','k')
area([0:.04:.98],real_counts,'FaceColor','k','FaceAlpha',.75,'LineWidth',1.5,'EdgeColor','k')
%histogram(sum(same_cohort,2)/100,0:.02:1,'FaceColor','k','Normalization','probability')
c = colorbar; colormap(map)

c.Ticks = [0:.2:1]
c.TickLabels = sim_thresh
c.Label.String = 'Similarity Threshold'
set(gca,'FontSize',18)
set(gcf,'Color','w')
xlabel('Mean Probability')
ylabel('Frequency')