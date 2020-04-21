function neuron = classify_nodes(neuron)
    node_sums = sum(neuron.skeleton_data.adj_uw,2);
    node_sums(node_sums == 1) = 0;
    node_sums(node_sums == 2) = 1;
    node_sums(node_sums > 2) = 2;
    neuron.skeleton_data.Node_Class = node_sums;
end