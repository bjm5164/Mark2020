%% Calculate synapse densities
% Calculate densities for a range of values
load mark2020_neurons_temporal_cohorts.mat
sigmas = [1000:1000:5000]
for jj = 1:length(sigmas)
    for i = 1:length(nl)
        input_densities = synapse_density(nl(i),2,sigmas(jj));
        output_densities = synapse_density(nl(i),1,sigmas(jj));
        nl(i).output_density(:,jj) = output_densities;
        nl(i).input_density(:,jj) = input_densities; clear input_densities and output_densities
    end
end

%% 
for kk = 2:5
    for i = 1:length(nl)
        for j = 1:length(nl)
            [~,n_real,n_possible] = filling_fraction(nl(i),nl(j),kk);
            real_adjacency{kk}(i,j) = n_real; clear n_real
            potential_adjacency{kk}(i,j) = n_possible; clear n_possible
        end
    end
end
%%
adj = get_adjacency(nl,0)
binary_potential = potential_adjacency
binary_adj = adj

binary_potential(binary_potential>1) = 1
binary_adj(binary_adj>1) = 1

corr2(binary_potential,binary_adj)


%% histogram of filling fractions

ff = adj./potential_adjacency
ff(isnan(ff)) = -1

histogram(ff(:))

idx = find(ff>1)
[ii,jj] = ind2sub(size(adj),idx)
figure; hold on
plot_neurons(nl(18),'m',1,3,1,1,0)
plot_neurons(nl(22),'g',1,3,1,1,0)