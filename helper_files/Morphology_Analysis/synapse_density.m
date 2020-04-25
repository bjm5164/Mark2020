function densities = synapse_density(neuron,direction,sigma)

    g_ud = graph((neuron.skeleton_data.Adj_dir+transpose(neuron.skeleton_data.Adj_dir)))
    input_nodes = arrayfun(@(x) neuron.skeleton_data.Tree2Ind(neuron.Inputs.treenodeID(x)),1:length(neuron.Inputs.treenodeID));
    distance_mat = distances(g_ud);
    sigma = 5000
    tic
    for i = 1:length(neuron.skeleton_data.Adj_dir)
        for j = 1:length(input_nodes)
            d_ij_in(i,j) = distance_mat(i,input_nodes(j));
            kernel_in(i,j) = exp((-1*(d_ij_in(i,j)^2)/(2*sigma^2)));
        end
        input_density(i) = sum(kernel_in(i,:));


    end
    toc

    density_range = [0:1:10,10:15:30]
    input_density_discrete = discretize(input_density,density_range)
    map = plasma(length(density_range))
    figure
    pg = plot_neurons(neuron,'k',1,3,1,1)
    for i = 1:length(input_density_discrete)
        highlight(pg,i,'NodeColor',map(input_density_discrete(i),:))
        highlight(pg,i,successors(neuron.skeleton_data.Graph,i),'EdgeColor',map(input_density_discrete(i),:))
    end
    axis equal
end


