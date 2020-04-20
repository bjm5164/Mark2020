function [neuron] = split_axon_dendrite(neuron)

if isempty(neuron.Outputs.treenodeID) == 0 & isempty(neuron.Inputs.treenodeID) == 0
    % Classify nodes into branches/ends/neither and get a binary adjacency
    % matrix for the skeleton.
    neuron = classify_nodes(neuron)
    adj = neuron.skeleton_data.Adj_dir;
    adj(adj>0) = 1;
    
    
    post_nodes = arrayfun(@(x) neuron.skeleton_data.tree2ind(neuron.Inputs.treenodeID(x)),1:length(neuron.Inputs.treenodeID));
    pre_nodes_c = arrayfun(@(x) repmat(neuron.skeleton_data.tree2ind(neuron.Outputs.treenodeID(x)),neuron.Outputs.polyadics(x),1),1:length(neuron.Outputs.treenodeID),'UniformOutput',false);
    pre_nodes = cat(1,pre_nodes_c{:}); clear output_nodes_c
    total_post = length(post_nodes);

    branches = neuron.skeleton_data.Node_Idxs(neuron.skeleton_data.Node_Class == 2)
    nodes_of_interest = [pre_nodes;post_nodes';branches]

    %Flip the towards soma graph to an away from soma graph.  Now only
    %distal synapses will have non-infinite distances.  Then find the
    %shortest path between pre or post synaptic points and the nodes of
    %interest (branches/synapses).  Anything that is non-zero/infinite is
    %distal node i.
    towards = digraph(adj)
    away = flipedge(towards)
    for i = 1:length(nodes_of_interest)
        for j = 1:length(pre_nodes)
            [~,D,~] = shortestpath(away,nodes_of_interest(i),pre_nodes(j));
            D_Pre(i,j) = D;
        end
    end
    D_Pre(D_Pre == inf) = 0; clear D  
    D_Pre(D_Pre > 0) = 1;

    for i = 1:length(nodes_of_interest)
        for j = 1:length(post_nodes)
            [~,D,~] = shortestpath(away,nodes_of_interest(i),post_nodes(j));
            D_Post(i,j) = D; 
        end
    end
    D_Post(D_Post == inf) = 0; clear D
    D_Post(D_Post>0) = 1;

    % Sum the number of distal synapses
    Distal_Pre = sum(D_Pre,2)
    Distal_Post = sum(D_Post,2)
    % 
    centrifugal = arrayfun(@(x) (total_post - Distal_Post(x)) * Distal_Pre(x),1:length(nodes_of_interest));
    centripital = arrayfun(@(x) Distal_Post(x) * (total_post - Distal_Pre(x)),1:length(nodes_of_interest));

    sum_flow = centrifugal + centripital;

    max_sfc = nodes_of_interest(find(sum_flow == max(sum_flow)))
    max_sfc = unique(max_sfc)

    if numel(max_sfc) > 1
        dtr = distances(graph(neuron.skeleton_data.Adj),neuron.skeleton_data.tree2ind(neuron.skeleton_data.Soma_ID),max_sfc);
        max_sfc = max_sfc(find(min(dtr)));
    else
       
    end
    
    % flip the towards soma graph, and get all distal points.  
    frag_a = dfsearch(away,max_sfc);
    frag_b = setdiff([1:length(neuron.skeleton_data.Adj_dir)],frag_a);
    
    frag_a_io_fraction = sum(ismember(frag_a,pre_nodes))/sum(ismember(frag_a,post_nodes));
    frag_b_io_fraction = sum(ismember(frag_b,pre_nodes))/sum(ismember(frag_b,post_nodes));
    
    if frag_a_io_fraction > frag_b_io_fraction
        axon = frag_a;
        dendrite = frag_b;
    else
        axon = frag_b;
        dendrite = frag_a;
    end
elseif isempty(neuron.Outputs.treenodeID) == 1
    axon = []
    dendrite = neuron.skeleton_data.Node_Idxs
elseif isempty(neuron.Inputs.treenodeID) == 1
    dendrite = []
    axon = neuron.skeleton_data.Node_Idxs
else
end
    
neuron.skeleton_data.axon_dendrite_index = ismember([1:length(adj)],axon)

    figure
    pc = plot_neurons(neuron,'k',1,3,1,0)
    %highlight(pc,neuron.skeleton_data.Node_Idxs(neuron.skeleton_data.Node_Class == 2),'NodeColor','m','MarkerSize',5)
    highlight(pc,max_sfc,'NodeColor','m','MarkerSize',10)
    
   highlight(pc,axon,'NodeColor','r','MarkerSize',4)
   highlight(pc,dendrite,'NodeColor','b','MarkerSize',4)
   
end