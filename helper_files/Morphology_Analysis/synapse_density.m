function densities = synapse_density(neuron,direction,sigma)

    g_ud = graph((neuron.skeleton_data.Adj_dir+transpose(neuron.skeleton_data.Adj_dir)))
    
    if direction == 1
            syn_nodes = arrayfun(@(x) neuron.skeleton_data.Tree2Ind(neuron.Outputs.treenodeID(x)),1:length(neuron.Outputs.treenodeID));
    elseif direction == 2
            syn_nodes = arrayfun(@(x) neuron.skeleton_data.Tree2Ind(neuron.Inputs.treenodeID(x)),1:length(neuron.Inputs.treenodeID));

    elseif direction == 3
            pre_nodes = arrayfun(@(x) neuron.skeleton_data.Tree2Ind(neuron.Outputs.treenodeID(x)),1:length(neuron.Outputs.treenodeID));
            post_nodes = arrayfun(@(x) neuron.skeleton_data.Tree2Ind(neuron.Inputs.treenodeID(x)),1:length(neuron.Inputs.treenodeID));
            syn_nodes = [pre_nodes,post_nodes]
        
    else error('wrong direction')
    end
    distance_mat = distances(g_ud);

  
    for i = 1:length(neuron.skeleton_data.Adj_dir)
        for j = 1:length(syn_nodes)
            d_ij_in(i,j) = distance_mat(i,syn_nodes(j));
            kernel_in(i,j) = exp((-1*(d_ij_in(i,j)^2)/(2*sigma^2)));
        end
        densities(i) = sum(kernel_in(i,:));


    end
   

%     density_range = [0,1.22.^[1:20]]
%     input_density_discrete = discretize(densities,density_range)
%     map = plasma(length(density_range))
%     figure
%     pg = plot_neurons(neuron,'k',1,3,1,1,1,map(input_density_discrete(neuron.skeleton_data.Soma),:))
%     for i = 1:length(input_density_discrete)
% %          if densities(i) > 5          
% %            highlight(pg,i,'NodeColor','b','MarkerSize',2)
% %          else
% %            highlight(pg,i,'NodeColor','k','MarkerSize',2)
% %          end
%          highlight(pg,i,'NodeColor',map(input_density_discrete(i),:))
%          highlight(pg,i,successors(neuron.skeleton_data.Graph,i),'EdgeColor',map(input_density_discrete(i),:))
%     end
%     axis equal
end


